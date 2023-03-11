#' Perform linear combination of the polygenic risk scores
#'
#' This function perform a linear combination of the scores
#'
#' @param pheno_file Directory to the phenotype file
#' @param covariate_file Directory to file with covariate information (age,sex,PC1..10)
#' @param score_files_list A vector contain directories of the PGS to be read
#' @param trait_specific_score_file A filename contain PGS IDs of trait-specific to combine (PRSmix), one score per line
#' @param pheno_name Column name of the phenotype in pheno_file
#' @param isbinary True if this is binary
#' @param out Prefix of output
#' @param metascore Meta-information from PGS Catalog contain PGS id and trait names. Must contains information for ALL the scores (DEFAULT = NULL)
#' @param liabilityR2 TRUE if liability R2 should be reported (DEFAULT = FALSE)
#' @param IID_pheno Column name of IID of phenotype file (e.g IID, person_id)
#' @param covar_list A vector of of covariates, must exists as columns in covariate_file (DEFAULT = age, sex, PC1..10))
#' @param ncores Number of CPU cores for parallel processing (DEFAULT = 1)
#' @param is_extract_adjSNPeff TRUE if extract adjSNPeff, FALSE if only calculate the combined PRS. May consume extended memory (DEFAULT = FALSE)
#' @param original_beta_files_list The vector of directories to SNP effect sizes used to compute original PRSs (DEFAULT = FALSE)
#' @param train_size_list A vector of training sample sizes. If NULL, all 80% of the samples will be used (DEFAULT = NULL)
#' @param power_thres_list A vector of power thresholds to select scores (DEFAULT = 0.95)
#' @param pval_thres_list A vector of P-value thresholds to select scores (DEFAULT = 0.05)
#' @param read_pred_training TRUE if PRSs were assessed in the training set was already run and can be read from file (DEFAULT = FALSE)
#' @param read_pred_testing TRUE if PRSs were assessed in the testing set was already run and can be read from file (DEFAULT = FALSE)
#' @return This function will return several files including 
#' - The case counts (for binary trait), 
#' - The dataframe of training and testing sample split from the main dataframe, 
#' - The prediction accuracy for each PRS in the training and testing set, 
#' - The prediction accuracy assessed in the testing set of the best PRS selected from the training set,
#' - the AUC of the NULL model of only covariates, the best PGS, PRSmix and PRSmix+ (adjusted for covariates), 
#' - Odds Ratio of the best PGS, PRSmix and PRSmix+ (adjusted for covariates), 
#' - The mixing weights of the scores used in combination, 
#' - The adjusted SNP effects to estimate PRSmix and PRSmix+ (if is_extract_adjSNPeff=TRUE)
#' - Return 0 if no error
#'
#' @importFrom stats cor dnorm lm glm logLik pchisq qchisq qnorm as.formula coef nobs pnorm power predict sd
#' @importFrom data.table fread fwrite
#' @importFrom dplyr bind_rows select all_of mutate group_by summarise rowwise filter
#' @importFrom stringr str_to_title
#' @importFrom magrittr %>%
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom caret train trainControl
#' @importFrom utils head read.table
#' @export
combine_PRS = function(
	pheno_file,
	covariate_file,
	score_files_list,
	trait_specific_score_file,
	pheno_name,
	isbinary,
	out,
	metascore = NULL,
	liabilityR2 = F,
	IID_pheno = "IID",
	covar_list = c("age", "sex", paste0("PC", 1:10)),
	ncores = 1,
	is_extract_adjSNPeff = F,
	original_beta_files_list = NULL,
	train_size_list = NULL,
	power_thres_list = c(0.95),
	pval_thres_list = c(0.05),
	read_pred_training = FALSE,
	read_pred_testing = FALSE,
	debug = F
	) {

	options(datatable.fread.datatable=FALSE)

	writeLines("--- Reading covariate data ---")
	basic_data = fread(covariate_file)

	writeLines("--- Reading all polygenic risk scores ---")
	
	all_scores = NULL

	for (score_file_i in 1:length(score_files_list)) {

		print(score_file_i)
		score_file = score_files_list[score_file_i]
		dd = fread(score_file)
		idx = which(endsWith(colnames(dd), "_SUM"))[-1]

		dd_sub = dd[,c(2,idx)]
		print(dim(dd_sub))
		if (is.null(all_scores)) {
			all_scores = dd_sub
		} else
			all_scores = merge(all_scores, dd_sub, by="IID")
	}

	colnames(all_scores)[2:ncol(all_scores)] = substring(colnames(all_scores)[2:ncol(all_scores)], 1, nchar(colnames(all_scores)[2:ncol(all_scores)])-4)

	pgs_list = fread(trait_specific_score_file, header=F)[,1]
	pgs_list = intersect(pgs_list, colnames(all_scores))

	pheno = fread(pheno_file)
	idx = which(colnames(pheno) %in% c(IID_pheno, pheno_name))
	pheno = pheno[,idx]
	colnames(pheno) = c("IID", "trait")

	writeLines("--- Merging Phenotype and PRS files ---")

	pheno_prs = merge(pheno, all_scores, by="IID")
	pheno_prs_cov = merge(pheno_prs, basic_data, by.x="IID", by.y=IID_pheno)
	pheno_prs_cov = pheno_prs_cov[which(!is.na(pheno_prs_cov$trait)),]

	#######################
	
	irnt = function(x) return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
	rr = function(x,d=3) round(x,d)

	#######################

	null_train_size_list = F
	if (is.null(train_size_list)) {
		writeLines("--- Selecting 80% as training data ---")
		train_size_list = floor(0.8*nrow(pheno_prs_cov))
		null_train_size_list = T
	}

	out_save = out

	for (train_size in train_size_list) {
		
		# train_size = train_size_list[1]
		
		writeLines(paste("--- Using ", train_size, " individuals for training sample ---"))
		
		if (!null_train_size_list) out = paste0(out_save, "_train.", train_size)
		
		set.seed(1)
		train_idx = sample(1:nrow(pheno_prs_cov), floor(0.8*nrow(pheno_prs_cov)))
		remaining_idx = c(1:nrow(pheno_prs_cov))[-train_idx]
		
		if (isbinary) fwrite(as.data.frame(table(pheno_prs_cov$trait)), paste0(out, "_case_counts.txt"), row.names=F, sep="\t", quote=F)
		
		train_df = pheno_prs_cov[train_idx,]
		test_df = pheno_prs_cov[-train_idx,]

		#### custom sensitivity train percentage
		set.seed(1)
		train_idx_sub = sample(1:nrow(train_df), train_size)
		train_df = train_df[train_idx_sub,]
		#############

		if (!isbinary) {
			train_df$trait = irnt(train_df$trait)
			test_df$trait = irnt(test_df$trait)
		}

		cc = apply(train_df[,covar_list], 2, function(x) length(unique(x)))
		covar_list1 = covar_list[-which(cc<=5)]
		for (i in covar_list1) train_df[i] = as.numeric(scale(train_df[i]))
		for (i in covar_list1) test_df[i] = as.numeric(scale(test_df[i]))

		fwrite(train_df[,c(1,2)], paste0(out, "_train_df.txt"), row.names=F, quote=F, sep="\t")
		fwrite(test_df[,c(1,2)], paste0(out, "_test_df.txt"), row.names=F, quote=F, sep="\t")

		################################ training #################################

		writeLines("--- Evaluating PRS in training set ---")
		
		training_file = paste0(out, "_train_allPRS.txt")
		read_pred_training_1 = (read_pred_training & file.exists(training_file))
		
		testing_file = paste0(out, "_test_allPRS.txt")
		read_pred_testing_1 = (read_pred_testing & file.exists(testing_file))
		
		if (!read_pred_training_1) {
			sumscore = apply(train_df[,3:ncol(train_df)], 2, sum)
			idx = which(sumscore==0)
			idx2 = which(names(idx)=="sex")
			if (length(idx2)>0) idx = idx[-idx2]
			if (length(idx)>0) train_df = train_df[,-match(names(idx), colnames(train_df))]

			pgs_list_all = colnames(train_df)
			# pgs_list_all = pgs_list_all[which(startsWith(pgs_list_all, "PGS"))]
			pgs_list_all = pgs_list_all[which(pgs_list_all %in% colnames(all_scores)[2:ncol(all_scores)])]
			pred_acc_train_allPGS_summary = eval_multiple_PRS(train_df, pgs_list_all, covar_list, liabilityR2, alpha=0.05, isbinary=isbinary)

			# pred_acc_train_allPGS_summary1 = pred_acc_train_allPGS_summary
			# pred_acc_train_allPGS_summary = pred_acc_train_allPGS_summary1

			fwrite(pred_acc_train_allPGS_summary, training_file, row.names=F, sep="\t", quote=F)
		} else {
			writeLines(paste0("Reading training file: ", training_file))
			pred_acc_train_allPGS_summary = fread(training_file)
		}

		pred_acc_train_allPGS_summary = as.data.frame(pred_acc_train_allPGS_summary)

		pred_acc_train_trait_summary = pred_acc_train_allPGS_summary
		pred_acc_train_trait_summary = pred_acc_train_trait_summary[order(as.numeric(pred_acc_train_trait_summary$pval), decreasing=F),]
		head(pred_acc_train_trait_summary)

		pred_acc_train_allPGS_summary1 = pred_acc_train_allPGS_summary %>%
			filter(pgs %in% pgs_list)
		pred_acc_train_allPGS_summary1 = pred_acc_train_allPGS_summary1[order(as.numeric(pred_acc_train_allPGS_summary1$R2), decreasing=T),]
		bestPRS = pred_acc_train_allPGS_summary1[1,1]
		bestPRS_acc = eval_single_PRS(test_df, pheno = "trait", prs_name=bestPRS, covar_list=covar_list, liabilityR2, alpha=0.05, isbinary=isbinary)

		fwrite(bestPRS_acc, paste0(out, "_best_acc.txt"), row.names=F, sep="\t", quote=F)

		################################ testing #################################

		writeLines("--- Evaluating PRS in testing set ---")

		if (!read_pred_testing_1) {

			if (debug)	writeLines("Running NULL model ---- ")

			sumscore = apply(test_df[,3:ncol(test_df)], 2, sum)
			idx = which(sumscore==0)
			idx2 = which(names(idx)=="sex")
			if (length(idx2)>0) idx = idx[-idx2]
			if (length(idx)>0) test_df = test_df[,-match(names(idx), colnames(test_df))]

			pred_acc_test_trait = eval_multiple_PRS(test_df, pgs_list,  covar_list, liabilityR2, alpha=0.05, isbinary=isbinary)

			pred_acc_test_trait_summary = pred_acc_test_trait
			pred_acc_test_trait_summary = pred_acc_test_trait_summary[order(as.numeric(pred_acc_test_trait_summary$pval), decreasing=F),]
			head(pred_acc_test_trait_summary)

			fwrite(pred_acc_test_trait_summary, testing_file, row.names=F, sep="\t", quote=F)
		} else {
			writeLines(paste0("Reading testing file: ", testing_file))
			pred_acc_test_trait_summary = fread(testing_file)
		}

		###########################################################################
		
		if (isbinary) {
			x_train = as.matrix(train_df %>% select(all_of(c(bestPRS, covar_list)), -trait))
			y_train = as.vector(train_df$trait)
			train_data = data.frame(x_train,trait=y_train)

			x_test = as.matrix(test_df %>% select(all_of(c(bestPRS, covar_list)), -trait))
			y_test = as.vector(test_df$trait)
			test_data = data.frame(x_test,trait=y_test)

			train_tmp = train_data[,c("trait", c(bestPRS, covar_list))]
			train_tmp$trait = as.factor(train_tmp$trait)
			
			formula_null = as.formula(paste0("trait ~ ", paste0(covar_list, collapse="+")))
			formula_bestPRS = as.formula(paste0("trait ~ scale(", bestPRS, ")+", paste0(covar_list, collapse="+")))
			
			ctrl = trainControl(
				method = "repeatedcv",
				allowParallel = TRUE,
				number = 3,
				verboseIter = T)

			cl = makePSOCKcluster(ncores)
			registerDoParallel(cl)

			set.seed(123)
			model_prsmix_null = train(
				formula_null, data = train_tmp, method = "glmnet",
				trControl = ctrl, family = "binomial",
				tuneLength = 50, verbose=T
			)
			
			set.seed(123)
			model_prsmix_bestPRS = train(
				formula_bestPRS, data = train_tmp, method = "glmnet",
				trControl = ctrl, family = "binomial",
				tuneLength = 50, verbose=T
			)

			stopCluster(cl)

			test_data1 = test_data
			
			test_pred = predict(model_prsmix_null, test_data1, type = "prob")[,2]
			auc_ci = ci.auc(test_data1$trait, test_pred)
			auc_out = data.frame(method="null_model", auc=auc_ci[2], lowerCI=auc_ci[1], upperCI=auc_ci[3])
			fwrite(auc_out, paste0(out, "_auc_NULL.txt"), row.names=F, sep="\t", quote=F)
			
			test_pred = predict(model_prsmix_bestPRS, test_data1, type = "prob")[,2]
			auc_ci = ci.auc(test_data1$trait, test_pred)
			auc_out = data.frame(method="bestPGS", auc=auc_ci[2], lowerCI=auc_ci[1], upperCI=auc_ci[3])
			fwrite(auc_out, paste0(out, "_auc_BestPGS.txt"), row.names=F, sep="\t", quote=F)
			

			model1 = glm(formula_bestPRS, data=test_data1, family="binomial")
			model1s = summary(model1)
			mm = exp(model1s$coefficients[2,1])
			ll = exp(model1s$coefficients[2,1] - 1.97*model1s$coefficients[2,2])
			uu = exp(model1s$coefficients[2,1] + 1.97*model1s$coefficients[2,2])
			pval = format.pval(model1s$coefficients[2,4])
			writeLines(paste0("OR = ", rr(mm), " (", rr(ll), "-", rr(uu), "); P-value=", pval))			
			fwrite(data.frame(mm, ll, uu, pval), paste0(out, "_OR_BestPGS.txt"), row.names=F, sep="\t", quote=F)

		}

		pred_acc_train_trait_summary = pred_acc_train_allPGS_summary %>%
			filter(pgs %in% pgs_list)
		
		pred_acc_test_trait_summary_out = pred_acc_test_trait_summary

		writeLines("--- Iterating power and p-value parameters --- ")
		for (power_thres in power_thres_list)
			for (pval_thres in pval_thres_list) {

				# pval_thres = pval_thres_list[2]
				# power_thres = power_thres_list[1]
				
				writeLines(paste0("P = ", pval_thres))
				writeLines(paste0("Power = ", power_thres))
				writeLines("PRSmix:")

				topprs = pred_acc_train_trait_summary %>%
					filter(pval <= pval_thres & power >= power_thres)

				head(topprs)
				topprs = topprs$pgs
				print(length(topprs))
				
				start_time = Sys.time()
				if (length(topprs) == 0) {
					print("No high power trait-specific PRS for PRSmix")
				} else {

					if (!isbinary) {
						
						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						if (length(topprs) > 1)	sd_train = apply(as.data.frame(x_train[,topprs]), 2, sd, na.rm=T)
						x_train[,topprs] = scale(x_train[,topprs])					
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)

						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))

						train_tmp = train_data[,c("trait", topprs, covar_list)]

						if (length(topprs) == 1) {
							ww = ww_raw = c(1)
							names(ww) = names(ww_raw) = topprs
						} else {

							ctrl = trainControl(
								method = "repeatedcv",
								allowParallel = TRUE,
								number = 3,
								verboseIter = T)

							cl = makePSOCKcluster(ncores)
							registerDoParallel(cl)

							set.seed(123)
							model_prsmix = train(
								formula, data = train_tmp, method = "glmnet",
								trControl = ctrl,
								tuneLength = 50, verbose=T
							)

							stopCluster(cl)

							model_prsmix$bestTune
							ww = coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[,1][-1]
							ww_raw = ww
							ww = ww[which(!names(ww) %in% covar_list)]
							ww = ww / sd_train[match(names(ww), names(sd_train))]
							
							if (all(ww == 0)) {
								writeLines("No weight for PRS")
								ww = c(1)
								names(ww) = bestPRS_acc$pgs
								topprs = bestPRS_acc$pgs
							}
						}
						
						test_df1 = cbind(test_data, IID=test_df$IID)
						test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)

						res_lm1 = eval_single_PRS(test_df1, pheno = "trait", prs_name="newprs", covar_list=covar_list, liabilityR2 = liabilityR2, alpha=pval_thres, isbinary=isbinary)
						
						res_lm1$pgs = "PRSmix"
						res_lm1

					} else {
						
						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						if (length(topprs) > 1)	sd_train = apply(as.data.frame(x_train[,topprs]), 2, sd, na.rm=T)
						x_train[,topprs] = scale(x_train[,topprs])						
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)

						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))
						
						train_tmp = train_data[,c("trait", topprs, covar_list)]
						train_tmp$trait = as.factor(train_tmp$trait)

						if (length(topprs) == 1) {
							ww = ww_raw = c(1)
							names(ww) = names(ww_raw) = topprs
						} else {
							
							ctrl = trainControl(
								method = "repeatedcv",
								allowParallel = TRUE,
								number = 3,
								verboseIter = T)
							
							cl = makePSOCKcluster(ncores)
							registerDoParallel(cl)
							
							set.seed(123)
							model_prsmix = train(
								formula, data = train_tmp, method = "glmnet",
								trControl = ctrl, family = "binomial",
								tuneLength = 50, verbose=T
							)

							stopCluster(cl)

							model_prsmix$bestTune
							ww = coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[,1][-1]
							ww_raw = ww
							ww = ww[which(!names(ww) %in% covar_list)]
							ww = ww / sd_train[match(names(ww), names(sd_train))]

							if (all(ww == 0)) {
								writeLines("No weight for PRS")
								ww = c(1)
								names(ww) = bestPRS_acc$pgs
								topprs = bestPRS_acc$pgs
							}							
							
							test_data1 = test_data
							test_data1[,topprs] = scale(test_data[,topprs])
							
							test_pred = predict(model_prsmix, test_data1, type = "prob")[,2]
							auc_ci = ci.auc(test_data1$trait, test_pred)
							auc_out = data.frame(method="prsmix", auc=auc_ci[2], lowerCI=auc_ci[1], upperCI=auc_ci[3])
							fwrite(auc_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_auc_PRSmix.txt"), row.names=F, sep="\t", quote=F)
							
						}
						
						test_df1 = cbind(test_data, IID=test_df$IID)
						test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
						
						res_lm1 = eval_single_PRS(test_df1, pheno="trait", prs_name="newprs", covar_list=covar_list, liabilityR2 = liabilityR2, alpha=pval_thres, isbinary=isbinary)
						
						res_lm1$pgs = "PRSmix"
						res_lm1
						
						############## OR ###################
						
						ff = paste0("trait ~ scale(newprs) + ", paste0(covar_list, collapse="+"))
						model1 = glm(ff, data=test_df1, family="binomial")
						model1s = summary(model1)
						mm = exp(model1s$coefficients[2,1])
						ll = exp(model1s$coefficients[2,1] - 1.97*model1s$coefficients[2,2])
						uu = exp(model1s$coefficients[2,1] + 1.97*model1s$coefficients[2,2])
						pval = format.pval(model1s$coefficients[2,4])
						print(paste0("OR = ", rr(mm), " (", rr(ll), "-", rr(uu), "); P-value=", pval))
						fwrite(data.frame(mm, ll, uu, pval), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_OR_PRSmix.txt"), row.names=F, sep="\t", quote=F)

						####################################

					}

					fwrite(data.frame(c(topprs), ww), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmix.txt"), sep="\t", quote=F)
					
					if (is_extract_adjSNPeff) {
						mixing_weight_file = paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmix.txt")
						outfile = paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_adjSNPeff_PRSmix.txt")
						extract_adjSNPeff(mixing_weight_file, original_beta_files_list, outfile)
					}
					
					res_lm1_summary = res_lm1
					res_lm1_summary$pgs = "PRSmix"
					pred_acc_test_trait_summary_out = bind_rows(res_lm1, pred_acc_test_trait_summary)
					head(pred_acc_test_trait_summary_out)
					
					fwrite(pred_acc_test_trait_summary_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_test_summary_traitPRS_withPRSmix.txt"), row.names=F, sep="\t", quote=F)
					
					prs_out = test_df1 %>%
						select(IID, newprs)
					colnames(prs_out) = c("IID", "prsmix")
					
					fwrite(prs_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_prsmix.txt"), row.names=F, sep="\t", quote=F)
					
				}

				end_time = Sys.time()
				timerunning = difftime(end_time, start_time, units = "secs")[[1]]

				timedf = data.frame(pgs="PRSmix", npgs=length(ww_raw), time=timerunning)
				print(timedf)

				fwrite(timedf, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_time_PRSmix.txt"), row.names=F, sep="\t", quote=F)		

				############################

				writeLines("PRSmix+:")

				topprs = pred_acc_train_allPGS_summary %>%
					filter(pval <= pval_thres & power >= power_thres)
				topprs = topprs$pgs
				print(length(topprs))

				start_time = Sys.time()

				if (length(topprs) == 0) {
					print("No high power trait-specific PRS for PRSmix")
				} else {

					if (!isbinary) {

						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						if (length(topprs) > 1)	sd_train = apply(as.data.frame(x_train[,topprs]), 2, sd, na.rm=T)
						x_train[,topprs] = scale(x_train[,topprs])
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)

						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))

						train_tmp = train_data[,c("trait", topprs, covar_list)]

						if (length(topprs) == 1) {
							ww = c(1)
							names(ww) = topprs
						} else {

							ctrl = trainControl(
								method = "repeatedcv",
								allowParallel = TRUE,
								number = 3,
								verboseIter = T)

							cl = makePSOCKcluster(ncores)
							registerDoParallel(cl)

							set.seed(123)
							model_prsmix = train(
								formula, data = train_tmp, method = "glmnet",
								trControl = ctrl,
								tuneLength = 50, verbose=T
							)

							stopCluster(cl)

							model_prsmix$bestTune
							ww = coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[,1][-1]
							
							ww = ww[which(!names(ww) %in% covar_list)]
							ww_raw = ww
							ww = ww / sd_train[match(names(ww), names(sd_train))]
							
							if (all(ww == 0)) {
								writeLines("No weight for PRS")
								ww = c(1)
								names(ww) = bestPRS_acc$pgs
								topprs = bestPRS_acc$pgs
							}
						}

						test_df1 = cbind(test_data, IID=test_df$IID)
						test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
						res_lm = eval_single_PRS(test_df1, pheno="trait", prs_name="newprs", covar_list=covar_list, liabilityR2 = liabilityR2, alpha=pval_thres, isbinary=isbinary)
						res_lm$pgs = "PRSmix+"
						res_lm

						nonzero_w = names(ww[which(ww!=0)])

					} else {

						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						if (length(topprs) > 1)	sd_train = apply(as.data.frame(x_train[,topprs]), 2, sd, na.rm=T)
						x_train[,topprs] = scale(x_train[,topprs])
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)
						
						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))

						train_tmp = train_data[,c("trait", topprs, covar_list)]
						train_tmp$trait = as.factor(train_tmp$trait)

						if (length(topprs) == 1) {
							ww = c(1)
							names(ww) = topprs
						} else {
							
							ctrl = trainControl(
								method = "repeatedcv",
								allowParallel = TRUE,
								number = 3,
								verboseIter = T)
							
							cl = makePSOCKcluster(ncores)
							registerDoParallel(cl)
							
							set.seed(123)
							model_prsmix = train(
								formula, data = train_tmp, method = "glmnet",
								trControl = ctrl, family = "binomial",
								tuneLength = 50, verbose=T
							)
							
							stopCluster(cl)
							
							model_prsmix$bestTune
							ww = coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[,1][-1]
							ww = ww[which(!names(ww) %in% covar_list)]
							ww_raw = ww							
							
							ww = ww / sd_train[match(names(ww), names(sd_train))]
							
							if (all(ww==0)) {
								writeLines("No weight for PRS")
								ww = c(1)
								names(ww) = bestPRS_acc$pgs
								topprs = bestPRS_acc$pgs
							}
							
							test_data1 = test_data
							test_data1[,topprs] = scale(test_data[,topprs])
							
							test_pred = predict(model_prsmix, test_data1, type = "prob")[,2]
							auc_ci = ci.auc(test_data1$trait, test_pred)
							auc_out = data.frame(method="prsmixP", auc=auc_ci[2], lowerCI=auc_ci[1], upperCI=auc_ci[3])
							fwrite(auc_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_auc_PRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
						}

						test_df1 = cbind(test_data, IID=test_df$IID)
						test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
						res_lm = eval_single_PRS(test_df1, pheno="trait", "newprs", covar_list, liabilityR2 = liabilityR2, alpha=pval_thres, isbinary=isbinary)
						res_lm$pgs = "PRSmix+"
						res_lm
						
						nonzero_w = names(ww[which(ww!=0)])
						
						################### OR ######################
						ff = paste0("trait ~ scale(newprs) + ", paste0(covar_list, collapse="+"))
						model = glm(ff, data=test_df1, family="binomial")
						model1s = summary(model)
						mm = exp(model1s$coefficients[2,1])
						ll = exp(model1s$coefficients[2,1] - 1.97*model1s$coefficients[2,2])
						uu = exp(model1s$coefficients[2,1] + 1.97*model1s$coefficients[2,2])
						pval = format.pval(model1s$coefficients[2,4])
						print(paste0("OR = ", rr(mm), " (", rr(ll), "-", rr(uu), "); P-value=", pval))
						fwrite(data.frame(mm, ll, uu, pval), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_OR_PRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
						
					}

					pred_acc_test_trait_summary_out = bind_rows(res_lm, pred_acc_test_trait_summary_out)
					head(pred_acc_test_trait_summary_out)

					fwrite(pred_acc_test_trait_summary_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_test_summary_traitPRS_withPRSmixPlus.txt"), row.names=F, sep="\t", quote=F)

					prsmixplus = test_df1 %>% select(IID, newprs)

					fwrite(prsmixplus, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_prsmixPlus.txt"), row.names=F, sep="\t", quote=F)

					fwrite(data.frame(topprs, ww), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmixPlus.txt"), row.names=F, sep="\t", quote=F)					
					
					if (is_extract_adjSNPeff) {
						mixing_weight_file = paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmixPlus.txt")
						outfile = paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_adjSNPeff_PRSmixPlus.txt")
						extract_adjSNPeff(mixing_weight_file, original_beta_files_list, outfile)
					}
					
					fwrite(data.frame(topprs, ww_raw), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_raw_PRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
					
					##############################################
				}


				end_time = Sys.time()
				timerunning = difftime(end_time, start_time, units = "secs")[[1]]
				timedf = data.frame(pgs="PRSmix+", npgs=length(ww_raw), time=timerunning)
				print(timedf)

				fwrite(timedf, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_time_PRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
				

			}
	}
	writeLines("Finished")
	return(0)

}




#' Combine PGS
#'
#' This function combines multiple PGS into one score
#'
#' @param trait The name of the trait
#' @param anc Intended ancestry
#' @param pgslist PGS list of the trait
#' @param pheno_name Name of the phenotype column
#' @param isbinary True if this is binary
#' @param score_files_list A list contain PGS directory
#' @param basic_data_file Directory to file with covariate information (age,sex,PC1...,PC10)
#' @param metascore Meta-information contain PGS id and trait names
#' @param phenofile Directory to the phenotype file
#' @param score_pref Prefix of score files
#' @param IID_pheno Column name of IID of phenotype file (e.g IID, person_id)
#' @param out Output prefix
#' @param covar_list Array of covariates
#' @param ncores Number of cores to run
#' @param power_thres_list Power threshold to select scores
#' @param pval_thres P-value threshold to select scores
#' @param read_pred_training True if the training set PRS assessment already run and can be read from file (Default: FALSE)
#' @param read_pred_testing True if the training set PRS assessment already run and can be read from file (Default: FALSE)
#' @return Prediction accuracy of PRSmix
#' @export
combine_PGS = function(
	trait,
	anc,
	pgslist,
	pheno_name,
	isbinary,
	score_files_list,
	basic_data_file,
	metascore,
	phenofile,
	score_pref,
	out,
	IID_pheno = "person_id",
	covar_list = c("age", "sex", paste0("PC", 1:10)),
	ncores = 5,
	train_size_list = NULL,
	power_thres_list = c(0.95),
	pval_thres_list = c(0.05),
	read_pred_training = NULL,
	read_pred_testing = NULL
	) {

	options(datatable.fread.datatable=FALSE)


	writeLines("Read basic data")

	basic_data = fread(basic_data_file)

	writeLines("Read all pgs")

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

	pgs_list = fread(pgslist, header=F)[,1]
	pgs_list = intersect(pgs_list, colnames(all_scores))

	pheno = fread(phenofile)
	idx = which(colnames(pheno) %in% c(IID_pheno, pheno_name))
	pheno = pheno[,idx]
	colnames(pheno) = c("IID", "trait")

	writeLines("Merging files")


	pheno_prs = merge(pheno, all_scores, by="IID")
	pheno_prs_cov = merge(pheno_prs, basic_data, by.x="IID", by.y=IID_pheno)
	pheno_prs_cov = pheno_prs_cov[which(!is.na(pheno_prs_cov$trait)),]

	#######################

	irnt = function(x) return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
	
	#######################
	
	null_train_size_list = F
	if (is.null(train_size_list)) {
		train_size_list = floor(0.8*nrow(pheno_prs_cov))
		null_train_size_list = T
	}
	
	out_save = out
	
	for (train_size in train_size_list) {
		
		train_size = train_size_list[1]
		
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
		
		writeLines("Evaluate PRS in training set")
		
		if (is.null(read_pred_training) & file.exists(paste0(out, "_train_allPRS.txt"))) { read_pred_training_1 = T } else { read_pred_training_1 = F }
		if (is.null(read_pred_testing) & file.exists(paste0(out, "_test_summary_traitPRS.txt"))) {read_pred_testing_1 = T } else { read_pred_testing_1 = F }

		
		if (!read_pred_training_1) {
			sumscore = apply(train_df[,3:ncol(train_df)], 2, sum)
			idx = which(sumscore==0)
			idx2 = which(names(idx)=="sex")
			if (length(idx2)>0) idx = idx[-idx2]
			if (length(idx)>0) train_df = train_df[,-match(names(idx), colnames(train_df))]
			
			pgs_list_all = colnames(train_df)
			pgs_list_all = pgs_list_all[which(startsWith(pgs_list_all, "PGS"))]
			pred_acc_train_allPGS_summary = get_acc_prslist_optimized(train_df, pgs_list_all, covar_list, isbinary)
			
			fwrite(pred_acc_train_allPGS_summary, paste0(out, "_train_allPRS.txt"), row.names=F, sep="\t", quote=F)
		} else {
			pred_acc_train_allPGS_summary = fread(paste0(out, "_train_allPRS.txt"))
		}

		pred_acc_train_allPGS_summary = as.data.frame(pred_acc_train_allPGS_summary)

		pred_acc_train_trait_summary = pred_acc_train_allPGS_summary
		pred_acc_train_trait_summary = pred_acc_train_trait_summary[order(as.numeric(pred_acc_train_trait_summary$pval_partial_R2), decreasing=F),]
		head(pred_acc_train_trait_summary)
		
		pred_acc_train_allPGS_summary1 = pred_acc_train_allPGS_summary %>%
			filter(pgs %in% pgs_list)
		pred_acc_train_allPGS_summary1 = pred_acc_train_allPGS_summary1[order(as.numeric(pred_acc_train_allPGS_summary1$R2), decreasing=T),]
		bestPRS = pred_acc_train_allPGS_summary1[1,1]
		bestPRS_acc = eval_prs(test_df, bestPRS, covar_list, isbinary)
		fwrite(bestPRS_acc, paste0(out, "_best_acc.txt"), row.names=F, sep="\t", quote=F)

		################################ testing #################################

		writeLines("Evaluate PRS in testing set")

		if (!read_pred_testing_1) {
			pred_acc_test_trait = get_acc_prslist_optimized(test_df, pgs_list,  covar_list, isbinary)

			pred_acc_test_trait_summary = pred_acc_test_trait
			pred_acc_test_trait_summary = pred_acc_test_trait_summary[order(as.numeric(pred_acc_test_trait_summary$pval_partial_R2), decreasing=F),]
			head(pred_acc_test_trait_summary)


			fwrite(pred_acc_test_trait_summary, paste0(out, "_test_summary_traitPRS.txt"), row.names=F, sep="\t", quote=F)
		} else {
			pred_acc_test_trait_summary = fread(paste0(out, "_test_summary_traitPRS.txt"))
		}

		###########################################################################

		pred_acc_train_trait_summary = pred_acc_train_allPGS_summary %>%
			filter(pgs %in% pgs_list)
		
		pred_acc_test_trait_summary_out = pred_acc_test_trait_summary
		
		for (power_thres in power_thres_list)
			for (pval_thres in pval_thres_list) {

				pval_thres = pval_thres_list[1]
				power_thres = power_thres_list[1]
				
				writeLines("PRSmix:")

				topprs = pred_acc_train_trait_summary %>%
					filter(pval_partial_R2 <= pval_thres & power >= power_thres)

				head(topprs)
				topprs = topprs$pgs
				print(length(topprs))

				print(length(topprs))
				if (length(topprs) == 0) {
					print("No high power trait-specific PRS for PRSmix")
				} else {
					
					if (!isbinary) {
						
						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)
						
						# formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))
						
						train_tmp = train_data[,c("trait", topprs, covar_list)]
						# train_tmp$trait = as.factor(train_tmp$trait)
						
						if (length(topprs) == 1) {
							ww = 1;
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
							
						}
						test_df1 = test_df
						test_df1$newprs = as.matrix(test_df1[,c(topprs, covar_list)]) %*% as.vector(ww)
						
						res_lm1 = eval_prs(test_df1, "newprs", covar_list, isbinary)
						res_lm1$pgs = "PRSmix"
						res_lm1
						
					} else {	
						
						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)
						
						# formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))
						
						train_tmp = train_data[,c("trait", topprs, covar_list)]
						train_tmp$trait = as.factor(train_tmp$trait)
						
						if (length(topprs) == 1) {
							ww = 1;
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
						}
						test_df1 = test_df
						# test_df1 = train_df
						test_df1$newprs = as.matrix(test_df1[,c(topprs, covar_list)]) %*% as.vector(ww)
						
						res_lm1 = eval_prs(test_df1, "newprs", covar_list, isbinary)
						res_lm1$pgs = "PRSmix"
						res_lm1

						############## OR ###################
						ff = paste0("trait ~ scale(newprs) + ", paste0(covar_list, collapse="+"))
						model1 = glm(ff, data=test_df1, family="binomial")
						model1s = summary(model1)
						mm = exp(model1s$coefficients[2,1])
						ll = exp(model1s$coefficients[2,1] - 1.97*model1s$coefficients[2,2])
						uu = exp(model1s$coefficients[2,1] + 1.97*model1s$coefficients[2,2])
						print(paste0(mm, " (", ll, "-", uu, ")"))
						fwrite(data.frame(mm, ll, uu), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_OR_PRSmix.txt"), row.names=F, sep="\t", quote=F)
						
						####################################
						
					}

					fwrite(data.frame(c(topprs, covar_list), ww), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PGSmix.txt"), sep="\t", quote=F)
					
					res_lm1_summary = res_lm1
					res_lm1_summary$pgs = "PRSmix"
					pred_acc_test_trait_summary_out = bind_rows(res_lm1, pred_acc_test_trait_summary)
					head(pred_acc_test_trait_summary_out)

					fwrite(pred_acc_test_trait_summary_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_test_summary_traitPRS_withPRSmix.txt"), row.names=F, sep="\t", quote=F)

					prs_out = test_df1 %>% 
						select(IID, pred_acc_test_trait_summary_out[2,1], newprs)
					colnames(prs_out) = c("IID", "pgscat", "prsmix")
					
					fwrite(prs_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_prsmix.txt"), row.names=F, sep="\t", quote=F)
				
				}

				############################

				writeLines("PRSmix+:")
				
				topprs = pred_acc_train_allPGS_summary %>%
					filter(pval_partial_R2 <= pval_thres & power >= power_thres)
				# topprs = pred_acc_train_allPGS_summary %>% filter(pgs=="PGS001590")
				# topprs = pred_acc_train_allPGS_summary %>% filter(power >= power_thres)
				topprs = topprs$pgs
				print(length(topprs))
				
				
				if (length(topprs) == 0) {
					print("No high power trait-specific PRS for PRSmix")
				} else {
					
					if (!isbinary) {
						
						
						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)
						
						# formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))
						
						train_tmp = train_data[,c("trait", topprs, covar_list)]
						# train_tmp$trait = as.factor(train_tmp$trait)
						
						if (length(topprs) == 1) {
							ww = 1;
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
							
						}
						
						test_df1 = test_df
						test_df1$newprs = as.matrix(test_df1[,c(topprs, covar_list)]) %*% as.vector(ww)
						res_lm = eval_prs(test_df1, "newprs", covar_list, isbinary)
						res_lm$pgs = "PRSmix+"
						res_lm

						nonzero_w = names(ww[which(ww!=0)])
						
					} else {
						
						x_train = as.matrix(train_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_train = as.vector(train_df$trait)
						train_data = data.frame(x_train,trait=y_train)

						x_test = as.matrix(test_df %>% select(all_of(c(topprs, covar_list)), -trait))
						y_test = as.vector(test_df$trait)
						test_data = data.frame(x_test,trait=y_test)
						
						# formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
						formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+", paste0(covar_list, collapse="+")))
						
						train_tmp = train_data[,c("trait", topprs, covar_list)]
						train_tmp$trait = as.factor(train_tmp$trait)
						
						if (length(topprs) == 1) {
							ww = 1;
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
							
						}
						
						test_df1 = test_df
						test_df1$newprs = as.matrix(test_df1[,c(topprs, covar_list)]) %*% as.vector(ww)
						res_lm = eval_prs(test_df1, "newprs", covar_list, isbinary)
						res_lm$pgs = "PRSmix+"
						res_lm

						nonzero_w = names(ww[which(ww!=0)])						

						################### OR ######################
						ff = paste0("trait ~ scale(newprs) + ", paste0(covar_list, collapse="+"))
						model = glm(ff, data=test_df1, family="binomial")
						models = summary(model)
						mm = exp(models$coefficients[2,1])
						ll = exp(models$coefficients[2,1] - 1.97*models$coefficients[2,2])
						uu = exp(models$coefficients[2,1] + 1.97*models$coefficients[2,2])
						print(paste0(mm, " (", ll, "-", uu, ")"))
						fwrite(data.frame(mm, ll, uu), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_OR_PRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
						
						
					}
					
					
					pred_acc_test_trait_summary_out = bind_rows(res_lm, pred_acc_test_trait_summary_out)
					head(pred_acc_test_trait_summary_out)

					fwrite(pred_acc_test_trait_summary_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_test_summary_traitPRS_withPRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
					
					prsmixplus = test_df1 %>% select(IID, newprs)
					
					fwrite(prsmixplus, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_prsmixPlus.txt"), row.names=F, sep="\t", quote=F)
					
					
					fwrite(data.frame(c(topprs,covar_list), ww), paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PGSmixPlus.txt"), row.names=F, sep="\t", quote=F)
					
					pgs_annot = fread(metascore)
					pgs_annot_sig = pgs_annot %>% filter(`Polygenic Score (PGS) ID` %in% nonzero_w)
					pgs_annot_sig_df = pgs_annot_sig %>%
						select(`Polygenic Score (PGS) ID`, `Reported Trait`)

					pgs_annot_sig_df = pgs_annot_sig_df[order(pgs_annot_sig_df$`Reported Trait`),]

					reported_trait = data.frame(table(pgs_annot_sig_df$`Reported Trait`))

					reported_trait = reported_trait %>%
						rowwise() %>%
						mutate(Var1 = gsub("\\s*\\([^\\)]+\\)","",Var1)) %>%
						mutate(Var1 = gsub("\\s*\\[[^\\)]+\\]","",Var1)) %>%
						mutate(Var1 = str_to_title(Var1)) %>%
						mutate(Var1 = gsub("Hdl","HDL",Var1)) %>%
						mutate(Var1 = gsub("Ldl","LDL",Var1)) %>%
						mutate(Var1 = gsub("Bmi","BMI",Var1)) %>%
						group_by(Var1) %>%
						summarise(Freq = sum(Freq))
						
					
					reported_trait = reported_trait[order(reported_trait$Freq, decreasing=T),]

					fwrite(reported_trait, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_reportedTraits.txt"), row.names=F, sep="\t", quote=F)
					reported_trait
					dim(reported_trait)
					
					
					##############################################
					
					
				}
			}
	}
	writeLines("Finished")
	return(0)

}




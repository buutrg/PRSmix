#' Combine PGS
#'
#' This function combines multiple PGS into one score
#'
#' @param trait The name of the trait
#' @param anc Intended ancestry
#' @param pgslist PGS list of the trait
#' @param score_pref Prefix of score files
#' @param phenofile Directory to the phenotype file
#' @param basic_data_file Directory to covariate information
#' @param pheno_name Name of the phenotype column
#' @param isbinary True if this is binary
#' @param out Output prefix
#' @return Prediction accuracy of PRSmix
#' @export
combine_PGS = function(
	trait = "cad",
	anc = "eur",
	# pgslist = "~/data/allPGSid.txt00",
	pgslist = "cad_list.txt",
	pheno_name = "CAD_case",
	isbinary=T,
	score_files_list = "",
	basic_data_file = "/home/jupyter/data/phenotypes/aou_basic_data.csv",
	metascore = "~/data/pgs_all_metadata_scores.csv",
	phenofile = "/home/jupyter/data/phenotypes/CAD_revised.csv",
	score_pref = "AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_", 
	out = "test_cad_mixedPRS"
	) {

	options(datatable.fread.datatable=FALSE)


	writeLines("Read basic data")

	basic_data = fread(basic_data_file)

	writeLines("Read all pgs")


	all_scores = NULL

	# for (i in 1:length(sscore_file_list)) {
	for (score_file_i in 1:length(score_files_list)) {
		
		# score_file_i = 1
		print(score_file_i)
		score_file = score_files_list[score_file_i]
		# i=1
		# dd = fread(paste0("~/data/prs_all/", sscore_file_list[i]))
		# dd = fread(paste0("~/data/optimization/", trait, "/", sscore_file_list[i]))
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
	idx = which(colnames(pheno) %in% c("person_id", pheno_name))
	pheno = pheno[,idx]
	colnames(pheno) = c("IID", "trait")


	writeLines("Merging files")


	pheno_prs = merge(pheno, all_scores, by="IID")
	pheno_prs_cov = merge(pheno_prs, basic_data, by.x="IID", by.y="person_id")
	pheno_prs_cov = pheno_prs_cov[which(!is.na(pheno_prs_cov$trait)),]

	############################

	irnt = function(x) return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))

	minmaxscaler = function(x, na.rm = TRUE) {
		s = (x- min(x)) / (max(x)-min(x))
	  return(as.numeric(s))
	}

	#######################

	set.seed(1)
	train_idx = sample(1:nrow(pheno_prs_cov), floor(.8*nrow(pheno_prs_cov)))
	remaining_idx = c(1:nrow(pheno_prs_cov))[-train_idx]

	# set.seed(1)
	# valid_idx = sample(remaining_idx, floor(.5*length(remaining_idx)))
	# test_idx = remaining_idx[which(!(remaining_idx %in% valid_idx))]


	train_df = pheno_prs_cov[train_idx,]
	test_df = pheno_prs_cov[-train_idx,]
	
	set.seed(1)
	valid_idx = sample(1:nrow(train_df), floor(0.2*nrow(train_df)))
	valid_df = pheno_prs_cov[valid_idx,]

	if (!isbinary) {
		train_df$trait = irnt(train_df$trait)
		test_df$trait = irnt(test_df$trait)
		valid$trait = irnt(valid$trait)
	}


	cov_list = c("age", paste0("PC", 1:10))
	for (i in cov_list) train_df[i] = as.numeric(scale(train_df[i]))
	for (i in cov_list) test_df[i] = as.numeric(scale(test_df[i]))
	for (i in cov_list) valid_df[i] = as.numeric(scale(valid_df[i]))

	################################ training #################################

	pgs_list_all = colnames(train_df)
	pgs_list_all = pgs_list_all[which(startsWith(pgs_list_all, "PGS"))]
	pred_acc_train_allPGS_summary = get_acc_prslist_optimized(train_df, pgs_list_all, isbinary)
	
	# pred_acc_train_allPGS_summary = fread(paste0(out, "_", anc, "_train_allPRS.txt"))
	
	fwrite(pred_acc_train_allPGS_summary, paste0(out, "_", anc, "_train_allPRS.txt"), row.names=F, sep="\t", quote=F)

	pred_acc_train_trait_summary = pred_acc_train_allPGS_summary
	pred_acc_train_trait_summary = pred_acc_train_trait_summary[order(as.numeric(pred_acc_train_trait_summary$pval_partial_R2), decreasing=F),]
	head(pred_acc_train_trait_summary)


	################################ testing #################################

	pred_acc_test_trait = get_acc_prslist_optimized(test_df, pgs_list, isbinary)

	pred_acc_test_trait_summary = pred_acc_test_trait
	pred_acc_test_trait_summary = pred_acc_test_trait_summary[order(as.numeric(pred_acc_test_trait_summary$pval_partial_R2), decreasing=F),]
	head(pred_acc_test_trait_summary)

	fwrite(pred_acc_test_trait_summary, paste0(out, "_", anc, "_test_summary_traitPRS.txt"), row.names=F, sep="\t", quote=F)

	###########################################################################

	pred_acc_train_trait_summary = pred_acc_train_allPGS_summary %>%
		filter(pgs %in% pgs_list)
	
	
	### Linear regression: trait specific

	if (!isbinary) {
		
		topprs = pred_acc_train_trait_summary %>% filter(pval_partial_R2 < 0.05 / nrow(pred_acc_train_trait_summary))
		topprs = topprs$pgs
		
		x_train = as.matrix(train_df %>% select(all_of(topprs), -trait))
		y_train = as.vector(train_df$trait)
		train_data = data.frame(x_train,trait=y_train)

		x_test = as.matrix(test_df %>% select(all_of(topprs), -trait))
		y_test = as.vector(test_df$trait)
		test_data = data.frame(x_test,trait=y_test)
		
		if (length(topprs) == 0) {
			print("No trait-specific significance")
		} else {
			
			formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
			
			train_tmp = train_data[,c("trait", topprs)]
			# train_tmp = as.matrix(train_tmp)
			
			ctrl <- trainControl(method = "repeatedcv", number = 5, verboseIter = T)
			set.seed(123)
			model <- train(
			  formula, data = train_tmp, method = "glmnet", 
			  trControl = ctrl,
			  tuneLength = 50, verbose=T
			)
			model$bestTune
			coef(model$finalModel, model$bestTune$lambda)
			ww = coef(model$finalModel, model$bestTune$lambda)[,1][-1]
			
			test_df1 = test_df
			test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
			
			
			# null_res = eval_null(test_df1, isbinary)
			
			res_lm1 = eval_prs(test_df1, null_res_test, "newprs", isbinary)
			res_lm1$pgs = "PRSmix"
			res_lm1
			
			pred_acc_detail_all1 = data.frame(pred_acc_test_trait_detail, "PRSmix"=res_lm1$prec_acc$partial_R2)
			fwrite(pred_acc_detail_all1, paste0(out, "_test_detailed_traitPRS_withPRSmix.txt"), row.names=F, sep="\t", quote=F)
			
		}
		
	} else {	
		
		# topprs = pred_acc_train_trait_summary %>% filter(pnew < 0.05 / nrow(pred_acc_train_trait_summary))
		# topprs = pred_acc_train_trait_summary %>% filter(pval_partial_R2 < 0.05)
		# topprs = pred_acc_test_trait_summary %>% filter(pval_partial_R2 < 0.05)
		topprs = pred_acc_train_trait_summary %>% filter(power >= 0.95)
		topprs = topprs$pgs
		# topprs = pred_acc_train_trait_summary$pgs
		
		print(length(topprs))
		
		x_train = as.matrix(train_df %>% select(all_of(topprs), -trait))
		y_train = as.vector(train_df$trait)
		train_data = data.frame(x_train,trait=y_train)

		x_test = as.matrix(test_df %>% select(all_of(topprs), -trait))
		y_test = as.vector(test_df$trait)
		test_data = data.frame(x_test,trait=y_test)
		
		x_valid = as.matrix(valid_df %>% select(all_of(topprs), -trait))
		y_valid = as.vector(valid_df$trait)
		valid_data = data.frame(x_valid,trait=y_valid)
		
		formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
		
		# train_tmp = train_data[,c("trait", topprs)]
		train_tmp = train_data[,c("trait", topprs)]
		train_tmp$trait = as.factor(train_tmp$trait)
		# train_tmp = as.matrix(train_tmp)
		
		ctrl <- trainControl(method = "repeatedcv",
	                        number = 5,
	                        # repeats = 5,
	                        # savePredictions = TRUE,
	                        verboseIter = T,
	                        returnResamp = "all")
		
		set.seed(123)
		model <- train(
		  formula, data = train_tmp, method = "glmnet", 
		  trControl = ctrl, family = "binomial", 
		  # intercept=FALSE,
		  tuneLength = 50, verbose=T
		)
		model$bestTune
		coef(model$finalModel, model$bestTune$lambda)
		ww = coef(model$finalModel, model$bestTune$lambda)[,1][-1]
		
		test_df1 = test_df
		test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
		
		res_lm1 = eval_prs(test_df1, "newprs", isbinary)
		res_lm1$pgs = "PRSmix"
		res_lm1
		
		
		
		############## OR ###################
		
		model1 = glm(trait ~ scale(newprs) + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC6 + PC7 + PC8 + PC9 + PC10, data=test_df1, family="binomial")
		model1s = summary(model1)
		mm = exp(model1s$coefficients[2,1])
		ll = exp(model1s$coefficients[2,1] - 1.97*model1s$coefficients[2,2])
		uu = exp(model1s$coefficients[2,1] + 1.97*model1s$coefficients[2,2])
		paste0(mm, " (", ll, "-", uu, ")")
		
		####################################
				
			
	# 		train_tmp1 = train_tmp
	# 		train_tmp1[,topprs] = scale(train_tmp1[,topprs])
			
	# 		formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
	# 		formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+"), "+age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC6 + PC7 + PC8 + PC9 + PC10"))
			
			
	# 		mod = glm(formula, train_data, family="binomial")
			
	# 		rf <- randomForest(formula, data=train_tmp1, ntree=50, do.trace=T, proximity=T)
	# 		print(rf)
			
	# 		test_data1 = test_data
	# 		test_data1[,topprs] = scale(test_data1[,topprs])
	# 		y_pred <- predict(rf, test_data1)
			
	# 		test_df1$rf = y_pred
	# 		res_lm1 = eval_prs(test_df1, null_res_test, "rf", isbinary)
	# 		res_lm1$pgs = "PRSmix_RF"
			
	# 		"age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC6", "PC7", "PC8", "PC9", "PC10"
	# 		feat = c(topprs, "age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC6", "PC7", "PC8", "PC9", "PC10")
	# 		feat = c(topprs)
	# 		# verbose = 1, print evaluation metric
	# 		dtrain <- xgb.DMatrix(data = x_train[,feat], label = y_train)
	# 		bst <- xgboost(data = dtrain, max.depth = 1000, eta = 1, nthread = 2, nrounds = 100, objective = "binary:logistic", verbose = 2)
			
	# 		y_pred <- predict(bst, as.matrix(x_test[,feat]))
	# 		prediction <- as.numeric(y_pred > 0.5)
	# 		table(prediction, y_test)
	# 		table(prediction)
	# 		print(head(prediction))
			
	# 		library(pROC)
	# 		auc(prediction, y_test)
	# 		err <- mean(as.numeric(y_pred > 0.5) != y_test)
	# print(paste("test-error=", err))
	}

	res_lm1_summary = res_lm1
	res_lm1_summary$pgs = "PRSmix"
	pred_acc_test_trait_summary_out = bind_rows(res_lm1, pred_acc_test_trait_summary)
	head(pred_acc_test_trait_summary_out)

	fwrite(pred_acc_test_trait_summary_out, paste0(out, "_", anc, "_test_summary_traitPRS_withPRSmix.txt"), row.names=F, sep="\t", quote=F)

	prs_out = test_df1 %>% 
		select(IID, pred_acc_test_trait_summary_out[2,1], newprs)
	colnames(prs_out) = c("IID", "pgscat", "prsmix")
	
	fwrite(prs_out, paste0(out, "_", anc, "_prsmix.txt"), row.names=F, sep="\t", quote=F)
	
	############################

	pred_acc_train_allPGS_summary = as.data.frame(pred_acc_train_allPGS_summary)

	all_sig = pred_acc_train_allPGS_summary %>% filter(power >= 0.95)
	all_sigpgs = all_sig$pgs
	print(length(all_sigpgs))

	pgs_list_sig = all_sigpgs
	print(length(pgs_list_sig))

	if (!isbinary) {
		
		idx = which(as.numeric(pred_acc_train$pval_partial_R2) < 0.05 / nrow(pred_acc_train))
		if (length(idx) == 0) idx = 1
		topprs = pred_acc_train[idx,1]
		
		
			
		x_train = as.matrix(train_df %>% select(all_of(topprs), -trait))
		y_train = as.vector(train_df$trait)
		train_data = data.frame(x_train,trait=y_train)

		x_test = as.matrix(test_df %>% select(all_of(topprs), -trait))
		y_test = as.vector(test_df$trait)
		test_data = data.frame(x_test,trait=y_test)
		
		
		
		train_tmp = train_data[,c("trait", topprs)]
		# train_tmp = as.matrix(train_tmp)
		
		formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
		
		ctrl <- trainControl(method = "repeatedcv", number = 5, verboseIter = T)
		set.seed(123)
		model <- train(
		  formula, data = train_tmp, method = "glmnet", 
		  trControl = ctrl,
		  tuneLength = 50, verbose=T
		)
		model$bestTune
		coef(model$finalModel, model$bestTune$lambda)
		ww = coef(model$finalModel, model$bestTune$lambda)[,1][-1]
		
		test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
		res_lm = eval_prs(test_df1, "newprs", isbinary)
		res_lm$pgs = "PRSmix+"
		
		res_lm_summary = res_lm
		res_lm_detail = res_lm$prec_acc
		
		
		pred_acc_detail_all1 = data.frame(pred_acc_test_trait_detail, "PRSmix+"=res_lm$prec_acc$partial_R2)
		fwrite(pred_acc_detail_all1, paste0(out, "_test_detailed_traitPRS_withPRSmixPlus_", anc, ".txt"), row.names=F, sep="\t", quote=F)		
		
		
	} else {
		
		topprs = pgs_list_sig
			
		formula = as.formula(paste0("trait ~ ", paste0(topprs, collapse="+")))
		
		
		
			
		x_train = as.matrix(train_df %>% select(all_of(topprs), -trait))
		y_train = as.vector(train_df$trait)
		train_data = data.frame(x_train,trait=y_train)

		x_test = as.matrix(test_df %>% select(all_of(topprs), -trait))
		y_test = as.vector(test_df$trait)
		test_data = data.frame(x_test,trait=y_test)
		
		
		train_tmp = train_data[,c("trait", topprs)]
		train_tmp$trait = as.factor(train_tmp$trait)
		# train_tmp = as.matrix(train_tmp)
		
		ctrl <- trainControl(method = "repeatedcv",
	                        number = 5,
	                        # savePredictions = TRUE,
	                        verboseIter = T,
	                        returnResamp = "all")
		
		set.seed(123)
		model <- train(
		  formula, data = train_tmp, method = "glmnet", 
		  trControl = ctrl, family = "binomial", 
		  # intercept=FALSE,
		  tuneLength = 50, verbose=T
		)
		model$bestTune
		coef(model$finalModel, model$bestTune$lambda)
		ww = coef(model$finalModel, model$bestTune$lambda)[,1][-1]
		nonzero_w = names(ww[which(ww!=0)])
		
		test_df1 = test_df
		test_df1$newprs = as.matrix(test_df1[,topprs]) %*% as.vector(ww)
		res_lm = eval_prs(test_df1, "newprs", isbinary)
		res_lm$pgs = "PRSmix+"
		
		res_lm_summary = res_lm
		# res_lm_detail = res_lm$prec_acc
		
		################### OR ######################
		
		model = glm(trait ~ scale(newprs) + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC6 + PC7 + PC8 + PC9 + PC10, data=test_df1, family="binomial")
		models = summary(model)
		mm = exp(models$coefficients[2,1])
		ll = exp(models$coefficients[2,1] - 1.97*models$coefficients[2,2])
		uu = exp(models$coefficients[2,1] + 1.97*models$coefficients[2,2])
		paste0(mm, " (", ll, "-", uu, ")")
		
		
	}
	
	
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

	fwrite(reported_trait, paste0(out, "_reportedTraits_", anc, ".txt"), row.names=F, sep="\t", quote=F)
	reported_trait
	dim(reported_trait)
	
	
	# pred_acc_test_trait_summary_out = rbind(res_lm_summary, pred_acc_test_trait_summary_out)
	# head(pred_acc_test_trait_summary_out)


	# fwrite(pred_acc_test_trait_summary_out, paste0(out, "_test_summary_traitPRS_withPRSmixPlus.txt"), row.names=F, sep="\t", quote=F)



	pred_acc_test_trait_summary_out = bind_rows(res_lm_summary, pred_acc_test_trait_summary_out)
	head(pred_acc_test_trait_summary_out)

	fwrite(pred_acc_test_trait_summary_out, paste0(out, "_", anc, "_test_summary_traitPRS_withPRSmixPlus.txt"), row.names=F, sep="\t", quote=F)
	
	prsmixplus = test_df1 %>% select(IID, newprs)
	
	fwrite(prsmixplus, paste0(out, "_", anc, "_prsmixPlus.txt"), row.names=F, sep="\t", quote=F)
	
	return(prsmixplus)

}




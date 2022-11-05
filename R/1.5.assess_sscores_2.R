get_acc_prslist_sub1 = function(data_df, pgs_list, null_res=NULL, isbinary=F) {
	
	# pgs_list = "PGS000329"
	# data_df = train_df
	
	if (is.null(null_res)) {	
		print("Fitting null model")
		null_res = eval_null(data_df, isbinary)
		# write.table(null_res, paste0("null_", ancestry, "_train_logLik_50rep.txt"), row.names=F, sep="\t", col.names=F, quote=F)
		# write.table(null_res, "null_test_logLik_50rep.txt", row.names=F, sep="\t", col.names=F, quote=F)
	}
		
	print("Fitting PRS")
	pred_acc_test = NULL
	pred_acc_test_detail = NULL
	system("mkdir prs_acc")
	
	for (prs_i in 1:length(pgs_list)) {
		if (file.exists(paste0("prs_acc/", pgs_list[prs_i], "_", ancestry, "_detail.txt"))) next()
		# prs_i = 2
		print(prs_i)
		pred_acc_test_tmp = eval_prs(data_df, null_res, pgs_list[prs_i], isbinary=isbinary)
		print(pred_acc_test_tmp$summary)
		pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp$summary)	
		pred_acc_test_detail = rbind(pred_acc_test_detail, pred_acc_test_tmp$prec_acc$partial_R2)
		
		fwrite(data.frame(pred_acc_test_tmp$summary), paste0("prs_acc/", pgs_list[prs_i], "_", ancestry, "_summary.txt"), row.names=F, sep="\t")
		fwrite(data.frame(pred_acc_test_tmp$prec_acc$partial_R2), paste0("prs_acc/", pgs_list[prs_i], "_", ancestry, "_detail.txt"), row.names=F, sep="\t")
		
	}

	pred_acc_test = pred_acc_test[order(pred_acc_test$partial_R2, decreasing=T),]
	
	rownames(pred_acc_test_detail) = pgs_list
	colnames(pred_acc_test_detail) = paste0("rep", 1:50)
	pred_acc_test_detail = t(pred_acc_test_detail)
	pred_acc_test_detail = as.data.frame(pred_acc_test_detail)


	
	return(list(pred_acc_test=pred_acc_test, pred_acc_test_detail=pred_acc_test_detail))
}


#' Assess PRS prediction accuracy
#'
#' This function estimates PRS prediction accuracy for a list of PRS
#'
#' @param trait Trait name
#' @param plink_prefix The genotype file in plink format (bed/bim/fam)
#' @param pgslist PGS list to evaluate
#' @param phenofile Phenotype file
#' @param pheno_name Phenotype name in the phenotype file
#' @param ancestry Ancestry
#' @param isbinary Is this trait binary?
#' @param out Output prefix
#' @param step1 is this running step 1 to get accuracy of the null model
#' @return The dataframe of prediction accuracy of all scores
#' @export
assess_score = function(
	trait,
	plink_prefix,
	pgslist,
	phenofile,
	ancestry,
	pheno_name,
	isbinary,
	out,
	step1 = F
	) {

	# opt = data.frame(	
	# 	trait="breast_cancer",
	# 	pgslist="~/data/allPGSid.txt00",
	# 	plink_prefix="AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_afr",
	# 	phenofile="/home/jupyter/data/phenotypes/cancer_phenotypes.csv",
	# 	ancestry="afr",
	# 	isbinary=T,
	# 	pheno_name="BreastCAFemale",
	# 	out="breast_cancer_afr_EvalAllPGS_part0"	
	# 	)


	# opt = data.frame(
	# 	trait = "cad",
	# 	ancestry = "eur",
	# 	pgslist = "~/data/allPGSid.txt00",
	# 	# pgslist = "cad_list.txt",
	# 	pheno_name = "CAD_case",
	# 	step1 = T,
	# 	plink_prefix = "AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_afr_",
	# 	isbinary=T,
	# 	phenofile = "/home/jupyter/data/phenotypes/CAD_revised.csv",
	# 	out = "test_cad_mixedPRS"
	# 	)


	# opt = data.frame(
		
	# 	trait = "ldl",
	# 	pgslist = "~/data/allPGSid.txt07",
	# 	plink_prefix = "AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_eur",
	# 	phenofile = "/home/jupyter/data/phenotypes/final_lipid_phenotypes.csv",
	# 	ancestry = "eur",
	# 	pheno_name = "LDL",
	# 	isbinary = F,
	# 	out = "ldl_eur_EvalAllPGS_part7"
		
	# 	)

	# opt = data.frame(
	# 	trait = "ovarian_cancer",
	# 	ancestry = "eur",
	# 	pgslist = "ovarian_cancer_list.txt",
	# 	pheno_name="OvarianCA",
	# 	isbinary=T,
	# 	phenofile = "~/data/phenotypes/cancer_phenotypes.csv",
	# 	out = "eval_ovarian_cancer_pgsmixpp"
	# 	)

	# opt = data.frame(
	# 	trait = "bmi",
	# 	ancestry = "eur",
	# 	pgslist = "bmi_list.txt",
	# 	pheno_name="value_as_number",
	# 	isbinary=F,
	# 	phenofile = "~/data/phenotypes/bmi.csv",
	# 	out = "eval_bmi_pgsmixpp"
	# 	)

	# trait = trait
	# pgslist = pgslist
	# isbinary = isbinary
	# ancestry = ancestry
	# phenofile = phenofile
	# pheno_name  = pheno_name
	# out = out

	print(opt)

	writeLines("Read basic data")

	basic_data = fread("/home/jupyter/data/phenotypes/aou_basic_data.csv")

	# sscore_file_list = list.files("~/data/prs_all/")
	sscore_file_list = list.files(paste0("~/data/optimization/", trait))
	sscore_file_list = sscore_file_list[which(endsWith(sscore_file_list, "sscore") & startsWith(sscore_file_list, plink_prefix))]

	writeLines("Read all pgs")


	all_scores = NULL

	for (i in 1:length(sscore_file_list)) {
		print(i)
		# i=1
		# dd = fread(paste0("~/data/prs_all/", sscore_file_list[i]))
		dd = fread(paste0("~/data/optimization/", trait, "/", sscore_file_list[i]))
		idx = which(endsWith(colnames(dd), "_SUM"))[-1]
		
		dd_sub = dd[,c(2,idx)]
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

	###########################

	irnt = function(x) return(qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))

	minmaxscaler = function(x, na.rm = TRUE) {
		s = (x- min(x)) / (max(x)-min(x))
		return(as.numeric(s))
	}


	#######################

	set.seed(1)
	train_idx = sample(1:nrow(pheno_prs_cov), floor(.8*nrow(pheno_prs_cov)))
	remaining_idx = c(1:nrow(pheno_prs_cov))[-train_idx]

	train_df = pheno_prs_cov[train_idx,]
	test_df = pheno_prs_cov[-train_idx,]


	if (!isbinary) {
		train_df$trait = irnt(train_df$trait)
		test_df$trait = irnt(test_df$trait)
	}

	cov_list = c("age", paste0("PC", 1:10))
	for (i in cov_list) train_df[i] = as.numeric(scale(train_df[i]))
	for (i in cov_list) test_df[i] = as.numeric(scale(test_df[i]))

	###########################

	print(paste0("null_", ancestry, "_train_logLik_50rep.txt"))

	if (step1) {
		null_res = eval_null(train_df, isbinary)
		write.table(null_res, paste0("null_", ancestry, "_train_logLik_50rep.txt"), row.names=F, sep="\t", col.names=F, quote=F)
		stop("Finish step 1")
	}

	null_res = NULL
	if (file.exists(paste0("null_", ancestry, "_train_logLik_50rep.txt"))) {
		null_res = fread(paste0("null_", ancestry, "_train_logLik_50rep.txt"))[,1]
	}

	# ################################ training #################################
	writeLines("Evaluating PRS")
	print(length(pgs_list))
	print(dim(train_df))
	print(dim(null_res))

	pred_acc_train_trait = get_acc_prslist_sub1(train_df, pgs_list, null_res, isbinary)

	pred_acc_train_trait_summary = pred_acc_train_trait$pred_acc_test
	pred_acc_train_trait_summary = pred_acc_train_trait_summary[order(as.numeric(pred_acc_train_trait_summary$pval_partial_R2), decreasing=F),]
	head(pred_acc_train_trait_summary)

	fwrite(pred_acc_train_trait_summary, paste0(out, "_", ancestry, "_train_summary_traitPRS.txt"), row.names=F, sep="\t", quote=F)


	pred_acc_train_trait_detail = pred_acc_train_trait$pred_acc_test_detail
	fwrite(pred_acc_train_trait_detail, paste0(out, "_", ancestry, "_train_detailed_traitPRS.txt"), row.names=F, sep="\t", quote=F)

	################################ testing #################################

	# pred_acc_test_trait = get_acc_prslist_sub1(test_df, pgs_list, null_res, isbinary)

	# pred_acc_test_trait_summary = pred_acc_test_trait$pred_acc_test
	# pred_acc_test_trait_summary = pred_acc_test_trait_summary[order(as.numeric(pred_acc_test_trait_summary$partial_R2), decreasing=F),]
	# head(pred_acc_test_trait_summary)


	# fwrite(pred_acc_test_trait_summary, paste0(out, "_test_summary_traitPRS.txt"), row.names=F, sep="\t", quote=F)



	# pred_acc_test_trait_detail = pred_acc_test_trait$pred_acc_test_detail
	# fwrite(pred_acc_test_trait_detail, paste0(out, "_test_detailed_traitPRS.txt"), row.names=F, sep="\t", quote=F)

	###########################################################################
	return(pred_acc_train_trait_summary)
}
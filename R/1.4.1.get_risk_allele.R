#' Get risk allele by consensus across multiple PRS snp effect panels
#'
#' This function get the risk-increasing allele based on consensus across PRS panels
#'
#' @param traitname Name of the trait
#' @param pgslist PGS list of the trait
#' @param pgsmetafile Metafile from PGS catalog (pgs_all_metadata_scores.csv)
#' @param freqfile Allele frequency file by plink2, ends with .afreq 
#' @param anc Ancestry
#' @param out Output
#' @return A dataframe of SNP and risk-increasing allele
#' @export
get_risk_allele = function(
	traitname = "cad",
	pgslist = "~/data/optimization/cad/cad_list.txt",
	pgsmetafile = "~/data/pgs_all_metadata_scores.csv",
	freqfile = "/home/jupyter/data/AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_eur.afreq",
	anc = "eur",
	out = "RiskAllele_cad_eur.txt"
	) {

	# opt = data.frame(
	# 	traitname = "cad",
	# 	pgslist = "~/data/optimization/cad/cad_list.txt",
	# 	plinkfile = "AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_eur",
	# 	anc = "eur",
	# 	out = "eval_cad_pgs.txt"
	# 	)

	# opt = data.frame(
	# 	traitname = "colorectal_cancer",
	# 	pgslist = "~/data/optimization/colorectal_cancer/colorectal_cancer_list.txt",
	# 	anc = "eur",
	# 	pgsmetafile = "~/data/pgs_all_metadata_scores.csv",
	# 	plinkfile = "AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_eur",
	# 	phenofile = "/home/jupyter/data/phenotypes/cancer_phenotypes.csv",
	# 	isbinary = T,
	# 	freqfile = "/home/jupyter/data/AoU_98K_WGS_QCed_callrate.0.9_hwe.1e-15_maf0.0001_eur.afreq",
	# 	out = "eval_colorectal_cancer_pgs.txt"
	# 	)

	# traitname = opt$traitname
	# pgslist = opt$pgslist
	# anc = opt$anc
	# plinkfile = opt$plinkfile
	# pgsmetafile = opt$pgsmetafile
	# freqfile = opt$freqfile
	# phenofile = opt$phenofile
	# isbinary = T
	# out = opt$out

	options(datatable.fread.datatable=FALSE)


	pgsinfo = data.table::fread(pgsmetafile)

	pgs_list = data.table::fread(pgslist, header=F)[,1]

	pgsinfo_trait = pgsinfo %>% filter(`Polygenic Score (PGS) ID` %in% pgs_list)

	head(pgsinfo_trait$`Polygenic Score (PGS) ID`)
	print(table(pgsinfo_trait$`PGS Publication (PGP) ID`))

	writeLines("reading freq file")
	freq = data.table::fread(freqfile)

	snp_weight_all = freq %>% select(ID, ALT, REF)
	snp_weight_all$BETA = 0
	colnames(snp_weight_all) = c("SNP", "A1", "A2", "BETA")

	tmp_list = rep(0, nrow(snp_weight_all))
	# freq1 = freq$ALT_FREQS[match(snp_weight_all$SNP, freq$ID)]
	# freq1 = sqrt(2*freq1*(1-freq))
	# snp_weight_all = NULL

	all_snps = NULL

	print(length(pgs_list))

	writeLines("reading snp effects")

	cc = 0
	for (prs_i in 1:length(pgs_list)) {
		
		# prs_i = 1
		print(prs_i)
		
		f = paste0("/home/jupyter/data/prs_aou/", pgs_list[prs_i], "_in_aou_", anc,".txt")
		if (file.exists(f) && file.size(f) > 0) {
			panel = data.table::fread(f)
			colnames(panel) = c("SNP", "A1", "BETA")
		} else {
			next()
		}
		
		cc = cc + 1
		# panel$BETA = scale(panel$BETA) * new_weight[prs_i]
		panel$BETA = panel$BETA
		
		if (is.null(snp_weight_all)) {
			snp_weight_all = panel %>% select(SNP, A1, BETA)
		} else {
			
			idx = match(panel$SNP, snp_weight_all$SNP)
			idx1 = which(is.na(idx))
			if (length(idx1) > 0) {
				print(paste0(length(idx1), " new snps"))
			# 	df_newsnp = panel[idx1,]
			# 	snp_weight_all = rbind(snp_weight_all, df_newsnp)
				panel = panel[-idx1, ]
			}
			
			idx = match(panel$SNP, snp_weight_all$SNP)
			same_a1 = panel$A1 == snp_weight_all$A1[idx]
			
			idx_notsame = which(!same_a1)
			
			idx2 = which(!same_a1)
			if (length(idx2)>0) {
				# print(length(idx2))
				panel[idx2,]$BETA = -panel[idx2,]$BETA
			}
			
			panel = panel[,c(1,3)]
			
			idx3 = match(panel$SNP, snp_weight_all$SNP)
			tmp_list1 = tmp_list
			tmp_list1[idx3] = panel$BETA
			
			
			all_snps = unique(c(all_snps, panel[,1]))
			
			snp_weight_all = cbind(snp_weight_all, tmp_list1)
			
			# snp_weight_all[,ncol(snp_weight_all)] = snp_weight_all[,ncol(snp_weight_all)] / freq1
			colnames(snp_weight_all)[ncol(snp_weight_all)] = pgs_list[prs_i]
			
			
		}
		
	}
	
	print(cc)
	snp_weight_all = snp_weight_all %>% select(-BETA)
	snp_weight_all = snp_weight_all[match(all_snps, snp_weight_all$SNP),]
	print(dim(snp_weight_all))

	snp_weight_all_save = snp_weight_all
	
	######################################################
	
	snp_weight_all = as.data.frame(snp_weight_all_save)
	
	cc = table(pgsinfo_trait$`PGS Publication (PGP) ID`)
	cc = which(cc > 1)
	print(length(cc))
	
	if (length(cc)>0) {
		for (i in 1:length(cc)) {
			print(i)
			studyid = names(cc)[i]
			idx = which(pgsinfo_trait$`PGS Publication (PGP) ID` == studyid)
			pgsid = pgsinfo_trait$`Polygenic Score (PGS) ID`[idx]
			idx1 = match(pgsid, colnames(snp_weight_all))
			subdf = snp_weight_all[,c(1,2,3,idx1)]
			
			subdf = apply(snp_weight_all[,idx1], 1, function(x) {
				if (all(x) == 0) return(0)
				if (any(x) > 0) return(1)
				if (any(x) < 0) return(-1)
				})
			
			snp_weight_all = snp_weight_all[,-idx1]
			snp_weight_all_tmp = data.frame(subdf)
			colnames(snp_weight_all_tmp) = studyid
			snp_weight_all = cbind(snp_weight_all, snp_weight_all_tmp)
			
		}
	}

	######################################################

	snp_weight_all_save = snp_weight_all

	######################################################

	snp_weight_all = snp_weight_all_save

	pos_sign = apply(snp_weight_all[,4:ncol(snp_weight_all)], 1, function(x) sum(x>0))
	neg_sign = apply(snp_weight_all[,4:ncol(snp_weight_all)], 1, function(x) sum(x<0))

	length(which(neg_sign>0))

	idx1 = which(neg_sign - pos_sign > 0)

	print(length(idx1))
	snp_weight_all$A1[idx1] = snp_weight_all$A2[idx1]
	snp_weight_all[idx1, 4:ncol(snp_weight_all)] = -snp_weight_all[idx1, 4:ncol(snp_weight_all)]

	snp_weight_all_out = snp_weight_all %>% select(SNP, A1)

	fwrite(snp_weight_all_out, out, row.names=F, sep="\t", quote=F, na=0)
	
	return(snp_weight_all_out)
}


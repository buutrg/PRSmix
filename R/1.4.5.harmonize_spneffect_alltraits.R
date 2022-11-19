#' Harmonize snp effects across PGS to the alternative allele
#'
#' This function harmonize SNP effect across PGS to the alternative allele
#'
#' @param freq_file Frequency file with ID, ALT, REF columns (e.g from PLINK2)
#' @param pgs_folder Directory to folder contain PGS
#' @param pgs_list File contain PGS id on each line
#' @param out Filename of the output
#' @return 0 if no error
#' @export
harmonize_snpeffect_toALT = function(
	freq_file, 
	pgs_folder, 
	pgs_list,
	out) {
	
	options(datatable.fread.datatable=FALSE)

	# opt = data.frame(
	# 	freq_file = "~/data/AoU_98K_WGS_QCed_callrate.0.95_NOhwe_maf0.001_eur.afreq",
	# 	pgs_folder = "~/data/prs_aou_eur",
	# 	pgs_list = "~/data/allPGSid.txt00",
	# 	out = "~/data/snp_weight_all_alltraits_callrate.0.95_NOhwe_maf0.001_eur_0.txt"
	# 	)



	# pgs_folder = opt$pgs_folder
	# freq_file = opt$freq_file
	# pgs_list = opt$pgs_list
	# out = opt$out


	pgs_list = fread(pgs_list, header=F)[,1]
	length(pgs_list)


	writeLines("reading freq file")
	freq = fread(freq_file)
	freq = freq %>% select(ID, ALT, REF)


	snp_weight_all = freq
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
	# for (prs_i in c(1,2,9)) {
		
		# prs_i = 1
		
		f = paste0(pgs_folder, "/", pgs_list[prs_i], ".txt")
		if (file.exists(f) && file.size(f) > 0) {
			panel = fread(f)
			colnames(panel) = c("SNP", "A1", "BETA")
		} else {
			next()
		}
		
		if (prs_i %% 50 == 0) print(prs_i)
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


	fwrite(snp_weight_all, out, row.names=F, sep=" ", quote=F, na=0)
	
	return(0)
}

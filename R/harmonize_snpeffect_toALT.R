#' Harmonize snp effects across PGS to the alternative allele
#'
#' This function harmonize SNP effect across PGS to the alternative allele
#'
#' @param ref_file Reference file contain SNP ID (ID), reference allele (REF) and alternative allele (ALT) columns (e.g allele frequency output --freq from PLINK2)
#' @param pgs_folder Directory to folder contain PGS
#' @param pgs_list File contain file names of single PGS on each line. The files must exist in the pgs_folder folder
#' @param out Filename of the output
#' @return A file is written containing SNP effects harmonized to the ALT allelle with columns SNP ID, ALT allele and BETAs
#' @export
harmonize_snpeffect_toALT = function(
	ref_file,
	snp_col = "SNP",
	a1_col = "A1",
	beta_col = "BETA",
	pgs_folder, 
	pgs_list,
	suffix = ".txt"
	out) {
	
	options(datatable.fread.datatable=FALSE)

	pgs_list = fread(pgs_list, header=F)[,1]
	length(pgs_list)

	writeLines("Reading freq file")
	freq = fread(ref_file)
	freq = freq %>% select(ID, ALT, REF)

	snp_weight_all = freq
	snp_weight_all$BETA = 0
	colnames(snp_weight_all) = c("SNP", "A1", "A2", "BETA")

	tmp_list = rep(0, nrow(snp_weight_all))
	all_snps = NULL
	print(length(pgs_list))

	writeLines("Reading snp effects")

	cc = 0
	for (prs_i in 1:length(pgs_list)) {
		
		f = paste0(pgs_folder, "/", pgs_list[prs_i], suffix)
		if (file.exists(f) && file.size(f) > 0) {
			panel = fread(f)
			panel = panel[,c(snp_col, a1_col, beta_col)]
			colnames(panel) = c("SNP", "A1", "BETA")
		} else {
			next()
		}
		
		if (prs_i %% 50 == 0) print(prs_i)
		cc = cc + 1
		panel$BETA = panel$BETA
		
		if (is.null(snp_weight_all)) {
			snp_weight_all = panel %>% select(SNP, A1, BETA)
		} else {
			
			idx = match(panel$SNP, snp_weight_all$SNP)
			idx1 = which(is.na(idx))
			if (length(idx1) > 0) {
				print(paste0(length(idx1), " new snps"))
				panel = panel[-idx1, ]
			}
			
			idx = match(panel$SNP, snp_weight_all$SNP)
			same_a1 = panel$A1 == snp_weight_all$A1[idx]
			
			idx_notsame = which(!same_a1)
			
			idx2 = which(!same_a1)
			if (length(idx2)>0) {
				panel[idx2,]$BETA = -panel[idx2,]$BETA
			}
			
			panel = panel[,c(1,3)]
			
			idx3 = match(panel$SNP, snp_weight_all$SNP)
			tmp_list1 = tmp_list
			tmp_list1[idx3] = panel$BETA
			
			
			all_snps = unique(c(all_snps, panel[,1]))			
			snp_weight_all = cbind(snp_weight_all, tmp_list1)			
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


merge_SNP_df = function(res_chunk) {
	res_chunk_fin = NULL
	for (subchunk_i in 1:length(res_chunk)) {

		cat(subchunk_i, "..")
		res_chunk_sub = res_chunk[[subchunk_i]]

		if (is.null(res_chunk_fin)) { res_chunk_fin = res_chunk_sub; next; }
		shared_snp = intersect(res_chunk_fin$SNP, res_chunk_sub$SNP)

		idx1 = match(shared_snp, res_chunk_fin$SNP)
		idx2 = match(shared_snp, res_chunk_sub$SNP)
		if (length(shared_snp)>0) res_chunk_fin[idx1,colnames(res_chunk_sub)] = res_chunk_sub[idx2, ]

		nonshared = which(!res_chunk_fin$SNP %in% shared_snp)
		if (length(nonshared)>0) res_chunk_fin = bind_rows(res_chunk_fin, res_chunk_sub[-idx2, ])
	}
	return(res_chunk_fin)
}


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
	suffix = ".txt",
	ncores = 1,
	chunk_size = 1,
	isheader = T,
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

	chunk_size = min(chunk_size, length(pgs_list))
	pgs_list_chunk_list = split(pgs_list, ceiling(seq_along(pgs_list)/chunk_size))

	res_chunk_all = NULL

	for (chunk_i in 1:length(pgs_list_chunk_list)) {

		# chunk_i = 1
		writeLines(paste0("Chunk ", chunk_i))

		pgs_list_chunk = pgs_list_chunk_list[[chunk_i]]

		res_chunk = mclapply(1:length(pgs_list_chunk), function(prs_i) {
			# prs_i = 2
			
			f = paste0(pgs_folder, "/", pgs_list[prs_i], suffix)

			if (file.exists(f) && file.size(f) > 0) {
				panel = fread(f, header = isheader)
				panel = panel[,c(snp_col, a1_col, beta_col)]
				colnames(panel) = c("SNP", "A1", "BETA")
			} else {
				return(NULL)
			}
			
			idx = match(panel$SNP, snp_weight_all$SNP)
			idx1 = which(is.na(idx))
			
			idx = match(panel$SNP, snp_weight_all$SNP)
			same_a1 = panel$A1 == snp_weight_all$A1[idx]
			
			idx_notsame = which(!same_a1)
			
			idx2 = which(!same_a1)
			if (length(idx2)>0) {
				panel[idx2,]$BETA = -panel[idx2,]$BETA
			}
			
			panel = panel[,c(1,3)]
			colnames(panel) = c("SNP", pgs_list_chunk[prs_i])
			return(panel)
		}, mc.cores = ncores)

		writeLines(paste0("Merging chunk ", chunk_i))
		res_chunk_fin = merge_SNP_df(res_chunk)
		if (is.null(res_chunk_all)) {
			res_chunk_all = res_chunk_fin
		} else {
			res_chunk_all = merge_SNP_df(list(res_chunk_all, res_chunk_fin))
		}
		print(dim(res_chunk_all))
	}

	res_chunk_all$A1 = snp_weight_all$A1[match(res_chunk_all$SNP, snp_weight_all$SNP)]
	res_chunk_all$A2 = snp_weight_all$A2[match(res_chunk_all$SNP, snp_weight_all$SNP)]
	res_chunk_all = res_chunk_all %>% 
		relocate(A1, .after=SNP) %>%
		relocate(A2, .after=A1)
	writeLines(paste0("Writing output file: ", out))
	fwrite(snp_weight_all, out, row.names=F, sep=" ", quote=F, na=0)
	
	return(0)
}



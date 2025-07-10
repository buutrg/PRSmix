
merge_SNP_df = function(res_chunk) {
	res_chunk_fin = NULL
	for (subchunk_i in 1:length(res_chunk)) {

		# subchunk_i = 2
		cat(subchunk_i, "..")
		res_chunk_sub = as.data.frame(res_chunk[[subchunk_i]])

		if (is.null(res_chunk_fin)) { res_chunk_fin = res_chunk_sub; next; }
		shared_snp = intersect(res_chunk_fin$SNP, res_chunk_sub$SNP)

		idx1 = match(shared_snp, res_chunk_fin$SNP)
		idx2 = match(shared_snp, res_chunk_sub$SNP)
		if (length(shared_snp)>0) res_chunk_fin[idx1,colnames(res_chunk_sub)] = res_chunk_sub[idx2, ]
		
		nonshared = which(!res_chunk_sub$SNP %in% shared_snp)
		if (length(nonshared)>0) {
			if (length(shared_snp)>0) {
				res_chunk_fin = bind_rows(res_chunk_fin, res_chunk_sub[-idx2, ])
			} else {
				res_chunk_fin = bind_rows(res_chunk_fin, res_chunk_sub)
			}
		}
	}
	return(res_chunk_fin)
}


merge_SNP_df_binary = function(res_chunk, left, right) {

	# print(paste0(left, "..", right))
	if (left == right) { return( as.data.frame(res_chunk[[left]])) }
	if (left < right) {
		mid = floor((left + right) / 2)
		left_df = as.data.frame(merge_SNP_df_binary(res_chunk, left, mid))
		right_df = as.data.frame(merge_SNP_df_binary(res_chunk, mid+1, right))
		
		# left = res_chunk[[10]]
		# right = res_chunk[[13]]

		shared_snp = intersect(left_df$SNP, right_df$SNP)

		idx1 = match(shared_snp, left_df$SNP)
		idx2 = match(shared_snp, right_df$SNP)
		if (length(shared_snp)>0) left_df[idx1,colnames(right_df)] = right_df[idx2, ]
		
		nonshared = which(!right_df$SNP %in% shared_snp)
		if (length(nonshared)>0) {
			if (length(shared_snp)>0) {
				left_df = bind_rows(left_df, right_df[-idx2, ])
			} else {
				left_df = bind_rows(left_df, right_df)
			}
		}
		return(as.data.frame(left_df))
	}

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

	pgs_list = fread(pgs_list, header=F)

	if (ncol(pgs_list)==2) {
		if (any(is.na(pgs_list[,2])) | any(duplicated(pgs_list[,2]))) 
			stop("There are 2 columns in PGS list file and second column contain NA or duplicated values")
		writeLines("There are 2 columns in PGS list and the second column will be assigned as names")
		name_pgs_list = pgs_list[,2]
		pgs_list = pgs_list[,1]
	} else {
		name_pgs_list = pgs_list[,1]
		pgs_list = pgs_list[,1]
	}
	length(pgs_list)

	writeLines("Reading reference file")
	freq = fread(ref_file, verbose=F)
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
	pgs_list_chunk_list_name = split(name_pgs_list, ceiling(seq_along(pgs_list)/chunk_size))

	res_chunk_all = NULL

	for (chunk_i in 1:length(pgs_list_chunk_list)) {

		# chunk_i = 2
		writeLines(paste0("Chunk ", chunk_i))

		pgs_list_chunk = pgs_list_chunk_list[[chunk_i]]
		pgs_list_chunk_name = pgs_list_chunk_list_name[[chunk_i]]

		res_chunk = mclapply(1:length(pgs_list_chunk), function(prs_i) {
			# prs_i = 28
			
			f = paste0(pgs_folder, "/", pgs_list_chunk[prs_i], suffix)

			if (file.exists(f) && file.size(f) > 0) {
				panel = fread(f, header = isheader, verbose=F)
				if (nrow(panel) == 0) return(NULL)
				panel = panel[,c(snp_col, a1_col, beta_col)]
				colnames(panel) = c("SNP", "A1", "BETA")
			} else {
				return(NULL)
			}
			
			idx = match(panel$SNP, snp_weight_all$SNP)			
			same_a1 = panel$A1 == snp_weight_all$A1[idx]
			
			idx_notsame = which(!same_a1)
			
			idx2 = which(!same_a1)
			if (length(idx2)>0) {
				panel[idx2,]$BETA = -panel[idx2,]$BETA
			}
			
			panel = panel[,c(1,3)]
			# colnames(panel) = c("SNP", pgs_list_chunk[prs_i])
			colnames(panel) = c("SNP", pgs_list_chunk_name[prs_i])
			return(panel)
		}, mc.cores = ncores)

		res_chunk = res_chunk[which(!sapply(res_chunk, is.null))]

		writeLines(paste0("Merging chunk ", chunk_i))

		if (length(res_chunk) == 0) next;
		
		# res_chunk1 = res_chunk[10:20]
		# res_chunk_fin = merge_SNP_df(res_chunk1)
		# res_chunk_fin = merge_SNP_df_binary(res_chunk1, 1, length(res_chunk1))		
		# res_chunk_fin = merge_SNP_df(res_chunk)
		res_chunk_fin = merge_SNP_df_binary(res_chunk, 1, length(res_chunk))
		
		if (is.null(res_chunk_all)) {
			res_chunk_all = res_chunk_fin
		} else {
			res_chunk_all = merge_SNP_df_binary(list(res_chunk_all, res_chunk_fin), 1, 2)
		}
		# if ("PGS000010" %in% colnames(res_chunk_all)) print(min(res_chunk_all$PGS000010, na.rm=T))

		print(dim(res_chunk_all))
	}

	res_chunk_all$A1 = snp_weight_all$A1[match(res_chunk_all$SNP, snp_weight_all$SNP)]
	res_chunk_all$A2 = snp_weight_all$A2[match(res_chunk_all$SNP, snp_weight_all$SNP)]
	res_chunk_all = res_chunk_all %>% 
		relocate(A1, .after=SNP) %>%
		relocate(A2, .after=A1)

	res_chunk_all[is.na(res_chunk_all)] = 0

	writeLines(paste0("Writing output file: ", out))
	fwrite(res_chunk_all, out, row.names=F, sep=" ", quote=F, na=0)
	
	return(0)
}



#' Estimate PRS
#'
#' This function estimate PRS with a reference scores
#'
#' @param geno Genotype file in plink format (bed/bim/fam).
#' @param weight_file The per-allele SNP effect
#' @param out Name of output file.
#' @param start_col Index of starting column of SNP effect sizes to estimate PRS in the weight file (DEFAULT = 4).
#' @return Execute plink. Return 0 if successful.
#' @export
compute_PRS = function(
	geno,
	weight_file,
	out,
	plink2_path = NULL,
	start_col = 4,
	ispgen = F) {
	
	n = length(read.table(pipe(paste0("head -n1 ", weight_file)), header=F)[1,])

	if (is.null(plink2_path)) {
		download.file("https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240105.zip", "plink.zip")
		unzip("plink.zip")
		plink2_path = "./plink2"
	}

	if (!ispgen) {
		cmd = system(paste0(plink2_path, " --bfile ", geno, " --score ", weight_file, " cols=+scoresums no-mean-imputation header-read --score-col-nums ", start_col, "-", n, " --out ", out))
	} else {
		cmd = system(paste0(plink2_path, " --pfile ", geno, " --score ", weight_file, " cols=+scoresums no-mean-imputation header-read --score-col-nums ", start_col, "-", n, " --out ", out))
	}
	
	print(cmd)
	system(cmd)
	
	return(0)
}

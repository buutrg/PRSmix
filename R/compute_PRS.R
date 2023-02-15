#' Estimate PRS
#'
#' This function estimate PRS with a reference scores
#'
#' @param plink_exe The weight file.
#' @param geno Genotype file in plink format (bed/bim/fam).
#' @param weight_file The per-allele SNP effect
#' @param out Name of output file.
#' @param start_col Index of starting column to estimate PRS in the reference file (DEFAULT = 4).
#' @return Execute plink. Return 0 if successful.
#' @export
compute_PRS = function(
	plink_exe,
	geno,
	weight_file,
	out,
	start_col = 4) {
	
	n = length(read.table(pipe(paste0("head -n1 ", weight_file)), header=F)[1,])
	
	download.file("https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_20221024.zip", "plink.zip")
	unzip("plink.zip")
	
	cmd = system(paste0(plink_exe, " --bfile ", geno, " --score ", weight_file, " cols=+scoresums no-mean-imputation header-read --score-col-nums ", start_col, "-", n, " --out ", out))
	
	print(cmd)
	system(cmd)
	
	return(0)
}

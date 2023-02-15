#' Estimate PRS
#'
#' This function estimate PRS with a reference scores
#'
#' @param wfile The weight file.
#' @param geno Genotype file in plink format (bed/bim/fam).
#' @param ref_file The reference risk-increasing allele file.
#' @param out Name of output file.
#' @param start_col Index of starting column to estimate PRS in the reference file (default = 4).
#' @return Execute plink. Return 0 if successful.
#' @export
compute_PRS = function(
	plink_exe,
	geno,
	ref_file,
	out,
	start_col = 4) {
	
	n = length(read.table(pipe(paste0("head -n1 ", ref_file)), header=F)[1,])
	
	cmd = system(paste0(plink_exe, " --bfile ", geno, " --score ", ref_file, " cols=+scoresums no-mean-imputation header-read --score-col-nums ", start_col, "-", n, " --out ", out))
	
	print(cmd)
	system(cmd)
	
	return(0)
}

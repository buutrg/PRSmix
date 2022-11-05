#' Harmonize snp effect with provided risk allele
#'
#' This function harmonize SNP effect according the panel of risk allele
#'
#' @param wfile The weight file
#' @param ref The reference risk-increasing allele file
#' @param out Name of output file
#' @return A dataframe of SNP and risk-increasing allele
#' @export
harmonize_snpeffect = function(
	wfile, 
	ref, 
	out) {

	# opt = data.frame(
	#   ref = "RiskAllele_breast_cancer_eur.txt",
	#   out = "snp_weight_allpgs_breast_cancer_eur",
	#   stringsAsFactor=F
	#   )

	# wfile = opt$wfile
	# out = opt$out
	ref = fread(ref)


	panel = fread(wfile)
	idx = match(ref$SNP, panel$SNP)

	notmatch_allele = which(ref$A1 != panel$A1[idx])

	idx1 = match(ref$SNP[notmatch_allele], panel$SNP)
	panel$A1[idx1] = ref$A1[notmatch_allele]
	panel[idx1, 4:ncol(panel)] = panel[idx1, 4:ncol(panel)] * -1

	fwrite(panel, out, row.names=F, sep=" ", quote=F, na=0)
	
	return(0)
	# return(panel)
}

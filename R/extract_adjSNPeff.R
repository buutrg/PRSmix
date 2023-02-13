#' Extract combined SNP effects
#'
#' This function extract the combined SNP effects
#'
#' @param mixing_weight_file The mixing weight files, 2 columns: 1) Score ID 2) Mixing weights
#' @param snp_eff_files_list A vector of names of file of original SNP weights
#' @param out Name of the output file
#' @return Adjusted SNP effects
#' @export
extract_adjSNPeff = function(
	mixing_weight_file,
    snp_eff_files_list,
    out
	) {

	options(datatable.fread.datatable=FALSE)
    
    mixing_weight = fread(mixing_weight_file)

    snp_eff = NULL
    for (snp_eff_files in 1:length(snp_eff_files_list)) {

        print(snp_eff_files)
        score_file = snp_eff_files_list[snp_eff_files]
        dd = fread(score_file, nrow=1)
        
        idx = which(colnames(dd) %in% mixing_weight[,1])
        if (length(idx)==0) next()
        
        idx = paste0("$", idx)
        cmd = paste0("awk '{print $1,$2,$3,", paste(idx, collapse=","), "}' ", score_file)
        dd = fread(cmd = cmd)
        
        idx = which(colnames(dd) %in% mixing_weight[,1])
        idx2 = match(colnames(dd)[idx], mixing_weight[,1])
        snp_eff_vec = as.matrix(dd[,idx]) %*% as.vector(mixing_weight[idx2,2])
        snp_eff_vec_tmp = data.frame(dd[,1:3], beta=snp_eff_vec)
        
        if (is.null(snp_eff)) {
            snp_eff = snp_eff_vec_tmp
        } else {
            idx_shared = intersect(snp_eff$SNP, snp_eff_vec_tmp$SNP)
            
            idx1 = match(idx_shared, snp_eff$SNP)
            idx2 = match(idx_shared, snp_eff_vec_tmp$SNP)
            snp_eff$beta[idx1] = snp_eff$beta[idx1] + snp_eff_vec_tmp$beta[idx2]
            snp_eff = rbind(snp_eff, snp_eff_vec_tmp[-idx2,])
        }
        
    }

    fwrite(snp_eff, out, row.names=F, sep="\t", quote=F)
    return(snp_eff)
}
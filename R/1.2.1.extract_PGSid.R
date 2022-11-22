#' Get PGS id
#'
#' This function extract PGS id from EFO code from PGS Catalog
#'
#' @param data directory to metadata score from PGS catalog (pgs_all_metadata_scores.csv)
#' @param trait The name of the trait
#' @param efo Ontology Trait ID (EFO from PGS catalog)
#' @return a dataframe of PGS id from PGS catalog
#' @export
extract_PGSid = function(
  data,
  trait,
  efo
  ) {

  options(datatable.fread.datatable=FALSE)

  efo_list = unlist(strsplit(opt$efo, split=","))

  pgsdata = fread(opt$data)

  pgsdata = pgsdata %>% 
    rowwise() %>%
    filter(sum(unlist(lapply(efo_list, function(x) grep(x, `Mapped Trait(s) (EFO ID)`))))>0)


  sum(unlist(lapply(efo_list, function(x) grep(x, pgsdata$`Mapped Trait(s) (EFO ID)`[1413]))))

  dim(pgsdata)
  
  res = pgsdata %>% select("Polygenic Score (PGS) ID")

  fwrite(res, paste0(opt$trait, "_list.txt"), row.names=F, col.names=F, quote=F)

  return(res)
}
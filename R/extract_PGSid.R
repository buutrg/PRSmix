#' Get PGS IDs from PGS Catalog using EFO id
#'
#' This function extract all PGS IDs (e.g. PGS000013, ...) from EFO code (e.g. EFO_0001645) from PGS Catalog
#'
#' @param data directory to metadata score from PGS catalog (pgs_all_metadata_scores.csv)
#' @param trait The name of the trait
#' @param efo Ontology Trait ID (EFO from PGS catalog)
#' @return a dataframe of PGS id from PGS catalog and a file with the PGS ID on each row
#' @export
extract_PGSid = function(
  data,
  trait,
  efo
  ) {

  options(datatable.fread.datatable=FALSE)

  efo_list = unlist(strsplit(efo, split=","))

  pgsdata = fread(data)

  pgsdata = pgsdata %>% 
    rowwise() %>%
    filter(sum(unlist(lapply(efo_list, function(x) grep(x, .data$`Mapped Trait(s) (EFO ID)`))))>0)

  sum(unlist(lapply(efo_list, function(x) grep(x, pgsdata$`Mapped Trait(s) (EFO ID)`[1413]))))

  dim(pgsdata)
  
  res = pgsdata %>% select("Polygenic Score (PGS) ID")

  fwrite(res, paste0(trait, "_list.txt"), row.names=F, col.names=F, quote=F)

  return(res)
}
#' Make dummy columns for categorical variables
#'
#' This function makes dummy columns for categorical variables
#'
#' @param df Dataframe with variables
#' @param cat_covar_list Categorical variables 
#' @return This function will generate a list containing a dataframe and its dummy columns for categorical variables 
#' - Return a list containing a dataframe and its dummy columns for categorical variables 
#'
#' @importFrom data.table fread fwrite
#' @importFrom dplyr bind_rows select all_of mutate group_by summarise rowwise filter
#' @importFrom fastDummies dummy_cols
#' @export
make_dummy_columns = function(df, cat_covar_list) {

    options(datatable.fread.datatable=FALSE)

    dummy_colnames = NULL

    for (cat_covar_i in 1:length(cat_covar_list)) {

        cat_covar = cat_covar_list[cat_covar_i]
        df = df[which(!is.na(df[,cat_covar])),]
        df = dummy_cols(df, select_columns = cat_covar, remove_selected_columns=T, remove_most_frequent_dummy=T)
        dummy_colnames_tmp = colnames(df)[which(startsWith(colnames(df), cat_covar))]
        dummy_colnames = c(dummy_colnames, dummy_colnames_tmp)
    }
    return(list(df=df, dummy_names=dummy_colnames))

}


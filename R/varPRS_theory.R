
######################## Get mixing weights with validation set

rr = function(x,digit=10) return(round(x,digit))



#' Evaluate PRS prediction accuracy
#'
#' This function get R2 for PRS accuracy
#'
#' @param data_df Data to assess prediction accuracy
#' @param prs_name PGS list of the trait
#' @param isbinary TRUE if binary and FALSE otherwise
#' @param covar_list Array of covariates
#' @param debug TRUE to verbose debugging
#' @return A dataframe for prediction accuracy of PRS and their power
#' @export
eval_prs = function(data_df, prs_name, covar_list, isbinary=F, debug=F) {
	
	prec_acc = NULL
	
	if (isbinary) {
		formula = as.formula(paste0("trait ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
		model_full = glm(formula, data=data_df, family="binomial")
		r_full = suppressWarnings(logLik(model_full, REML=FALSE))[1]
		
		formula = as.formula(paste0("trait ~ ", paste0(covar_list, collapse="+")))
		model_null = glm(formula, data=data_df, family="binomial")
		r_null = suppressWarnings(logLik(model_null, REML=FALSE))[1]
		
		m = r_full
		n = r_null
		N = nobs(model_full)
		cs = 1 - exp(-2/N * (m - n))
		nk = cs/(1 - exp(2/N * n))
		R2 = nk	
	} else {
		formula = as.formula(paste0("trait ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
		model_full = lm(formula, data=data_df)
		r_full = summary(model_full)$r.squared
		
		formula = as.formula(paste0("trait ~ ", paste0(covar_list, collapse="+")))
		model_null = lm(formula, data=data_df)
		r_null = summary(model_null)$r.squared
		
		N = nobs(model_full)
		R2 = r_full - r_null
	}
	
	alpha = 0.05
	NCP = N * R2 / (1-R2)
	power = 1-pnorm(qnorm(1-alpha/2)-NCP^0.5) + pnorm(qnorm(alpha/2)-NCP^0.5)
	
	vv = (4*R2*(1-R2)^2 *(N-2)^2) / ((N^2-1)*(N+3))

	se = sqrt(vv)
	lower_r2 = R2 - 1.97*se
	upper_r2 = R2 + 1.97*se
	pval = pchisq((R2/se)^2, df=1, lower.tail=F)
	
	r2_out = paste0(rr(R2), " (", rr(lower_r2), "-", rr(upper_r2), ")")
	
	return(data.frame(pgs=prs_name, R2=R2, partial_R2=r2_out, se=se, lowerCI=lower_r2, upperCI=upper_r2, pval_partial_R2=pval, power=power))
}


#' Get risk allele by consensus across multiple PRS snp effect panels
#'
#' This function get the risk-increasing allele based on consensus across PRS panels
#'
#' @param data_df Data to assess prediction accuracy
#' @param pgs_list PGS list of the trait
#' @param covar_list Array of covariates
#' @param isbinary TRUE if binary and FALSE otherwise
#' @return A dataframe for prediction accuracy of PRS and their power
#' @export
get_acc_prslist_optimized = function(data_df, pgs_list, covar_list, isbinary=F) {
	
	# data_df = test_df
	
	pred_acc_test = NULL
	for (prs_i in 1:length(pgs_list)) {
		
		if (prs_i %% 100 == 0) print(prs_i)
		
		prs_name = pgs_list[prs_i]
		pred_acc_test_tmp = eval_prs(data_df, prs_name, covar_list=covar_list, isbinary=isbinary)
		pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp)
	}
	
	pred_acc_test = pred_acc_test[order(pred_acc_test$R2, decreasing=T),]
	
	return(pred_acc_test)
}



######################## Get mixing weights with validation set

rr = function(x,digit=10) return(round(x,digit))

#' Evaluate prediction accuracy of a single PRS 
#'
#' This function get prediction accuracy for a single PRS
#'
#' @param data_df Data to assess prediction accuracy
#' @param prs_name PGS list of the trait
#' @param covar_list Array of covariates
#' @param isbinary TRUE if binary and FALSE otherwise (default = FALSE)
#' @param liabilityR2 TRUE if liability R2 should be reported (default = FALSE)
#' @param alpha Significance level to estimate power (default = 0.05)
#' @return A dataframe for prediction accuracy of a single PRS and their power with R2, R2 for output form, standard error, lower 95% CI, upper 95% CI, P-value and power
#' @export
eval_single_PRS = function(data_df, pheno = "trait", prs_name, covar_list, isbinary=F, liabilityR2=F, alpha=0.05) {
	
	if (isbinary & !liabilityR2) {
		data_df$trait = as.numeric(data_df$trait)
		formula = as.formula(paste0(pheno, " ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
		model_full = glm(formula, data=data_df, family="binomial")
		r_full = suppressWarnings(logLik(model_full, REML=FALSE))[1]
		
		formula = as.formula(paste0(pheno, " ~ ", paste0(covar_list, collapse="+")))
		model_null = glm(formula, data=data_df, family="binomial")
		r_null = suppressWarnings(logLik(model_null, REML=FALSE))[1]
		
		m = r_full
		n = r_null
		N = nobs(model_full)
		cs = 1 - exp(-2/N * (m - n))
		nk = cs/(1 - exp(2/N * n))
		R2 = nk
		
	} else {
		data_df$trait = as.numeric(data_df$trait)
		
		formula = as.formula(paste0(pheno, " ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
		model_full = lm(formula, data=data_df)
		r_full = summary(model_full)$r.squared
		
		formula = as.formula(paste0(pheno, " ~ ", paste0(covar_list, collapse="+")))
		model_null = lm(formula, data=data_df)
		r_null = summary(model_null)$r.squared
		
		N = nobs(model_full)
		R2 = r_full - r_null		
		
		if (liabilityR2) {
			N = nrow(data_df)
			K = mean(data_df$trait)
			R2 = R2 * K * (1-K) / (dnorm(qnorm(p=1-K, lower.tail=T))^2)
		}
		
	}
	
	NCP = N * R2 / (1-R2)
	power = 1-pnorm(qnorm(1-alpha/2)-NCP^0.5) + pnorm(qnorm(alpha/2)-NCP^0.5)
	
	vv = (4*R2*(1-R2)^2 *(N-2)^2) / ((N^2-1)*(N+3))
	
	se = sqrt(vv)
	lower_r2 = R2 - 1.97*se
	upper_r2 = R2 + 1.97*se
	pval = pchisq((R2/se)^2, df=1, lower.tail=F)
	
	r2_out = paste0(rr(R2,3), " (", rr(lower_r2,3), "-", rr(upper_r2,3), ")")
	
	return(data.frame(pgs=prs_name, R2=R2, se=se, lowerCI=lower_r2, upperCI=upper_r2, pval=pval, power=power))
}


#' Evaluate multiple PRSs from the given vector of PRS
#'
#' This function get prediction accuracy of the list of PRS
#'
#' @param data_df A dataframe where rows for samples, columns for PRS
#' @param pgs_list PGS list to evaluate, must be exist as columns in the dataframe
#' @param covar_list Array of covariates, must be exist as columns in the dataframe
#' @param liabilityR2 TRUE if liability R2 should be reported (default = FALSE)
#' @param alpha Significance level to estimate power (default = 0.05)
#' @param isbinary TRUE if binary and FALSE otherwise (default = FALSE)
#' @return A dataframe for prediction accuracy of multiple PRSss and their power with R2, R2 for output form, standard error, lower 95% CI, upper 95% CI, P-value and power
#'
#' @export
eval_multiple_PRS = function(data_df, pgs_list, covar_list, liabilityR2=F, alpha=0.05, isbinary=F) {
	
	writeLines("Eval all PGS")
	if (isbinary) print(table(data_df$trait))

	pred_acc_test = NULL
	for (prs_i in 1:length(pgs_list)) {
		
		if (prs_i %% 100 == 0) {
			writeLines(paste0("Evaluated ",  prs_i, " scores"))
		}
		
		prs_name = pgs_list[prs_i]
		pred_acc_test_tmp = eval_single_PRS(data_df=data_df, pheno="trait", prs_name=prs_name, covar_list=covar_list, isbinary=isbinary, liabilityR2=liabilityR2, alpha=alpha)
		pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp)
	}
	
	pred_acc_test = pred_acc_test[order(pred_acc_test$R2, decreasing=T),]
	
	return(pred_acc_test)
}


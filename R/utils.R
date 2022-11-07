
######################## Get mixing weights with validation set

rr = function(x,digit=10) return(round(x,digit))

getCI = function(x) {
	mm = mean(x)
	sdd = sd(x)
	n = length(x)
	# n = 1
	lowerCI = rr(mm - 1.97*sdd/sqrt(n))
	upperCI = rr(mm + 1.97*sdd/sqrt(n))
	return(paste0(rr(mm), " (", lowerCI, ";", upperCI, ")"))
}

wald = function(x) {
	mm = mean(x)
	sdd = sd(x)
	n = length(x)
	n = 1
	z = mm/(sdd/sqrt(n))
	p = 2*pnorm(q=z, lower.tail=FALSE) # two tail
	return(p)
}

eval_null = function(data_df, isbinary=F) {
	
	data_df = train_df
	prec_acc = NULL
	
	for (rep in 1:50) {
		if (rep %% 10 == 0) print(rep)
		set.seed(rep)
		data_df_sub = data_df[sample(1:nrow(data_df), floor(1*nrow(data_df)), replace=T),]
		# data_df_sub = data_df[sample(1:nrow(data_df), floor(0.9*nrow(data_df)), replace=T),]
		table(data_df_sub$trait)
		
		formula = as.formula(paste0("trait ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC6 + PC7 + PC8 + PC9 + PC10"))
		
		if (!isbinary) {
			model_null = lm(formula, data=data_df_sub)
			r_null = summary(model_null)$r.squared
		} else {
			model_null = glm(formula, data=data_df_sub, family="binomial")
			r_null = suppressWarnings(logLik(model_null, REML=FALSE))[1]
		}
		prec_acc = c(prec_acc, r_null)
	}
	
	return(prec_acc)
}

eval_prs = function(data_df, null_res, prs_name, isbinary=F, debug=F) {
	
	prec_acc = NULL
	for (rep in 1:50) {
		
		set.seed(rep)
		data_df_sub = data_df[sample(1:nrow(data_df), floor(1*nrow(data_df)), replace=T),]
		# data_df_sub = data_df[sample(1:nrow(data_df), floor(0.9*nrow(data_df)), replace=T),]
		table(data_df_sub$trait)
		
		formula = as.formula(paste0("trait ~ scale(", prs_name, ") + age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC6 + PC7 + PC8 + PC9 + PC10"))
		
		if (!isbinary) {
			model_full = lm(formula, data=data_df_sub)
			r_full = summary(model_full)$r.squared
			r_null = null_res[rep]
			partial_R2 = r_full - r_null
		} else {
			model_full = glm(formula, data=data_df_sub, family="binomial")
			
			m = suppressWarnings(logLik(model_full, REML=FALSE))[1]
			n = null_res[rep]
			N = nobs(model_full)
			cs = 1 - exp(-2/N * (m - n))
			nk = cs/(1 - exp(2/N * n))
			partial_R2 = nk
			# print(partial_R2)	
		}
		
		prec_acc = c(prec_acc, partial_R2)
	}
		
	prec_acc = data.frame(partial_R2=prec_acc)
	
	test_sig = apply(prec_acc, 2, wald)
	test_sig = as.numeric(test_sig)
	names(test_sig) = paste0("pval_", colnames(prec_acc))
	test_sig = t(data.frame(test_sig))
	
	prec_acc_ci = apply(prec_acc, 2, getCI)
	prec_acc1 = data.frame(pgs = prs_name, t(data.frame(prec_acc_ci)))
	prec_acc1 = data.frame(prec_acc1, test_sig)
	rownames(prec_acc1) = NULL
	return(list(summary=prec_acc1, prec_acc=prec_acc))
}



get_acc_prslist = function(data_df, pgs_list, null_res=NULL, isbinary=F) {
	
	# pgs_list = "PGS000329"
	# data_df = test_df
	
	if (is.null(null_res)) {	
		print("Fitting null model")
		null_res = eval_null(data_df, isbinary)
		# write.table(null_res_test, "null_train_logLik_50rep.txt", row.names=F, sep="\t", col.names=F, quote=F)
		# write.table(null_res_test, paste0("null_", anc,"test_logLik_50rep.txt"), row.names=F, sep="\t", col.names=F, quote=F)
	}
	
	print("Fitting PRS")
	pred_acc_test = NULL
	pred_acc_test_detail = NULL
	
	for (prs_i in 1:length(pgs_list)) {
		# prs_i = 2
		print(prs_i)
		pred_acc_test_tmp = eval_prs(data_df, null_res, pgs_list[prs_i], isbinary=isbinary)
		print(pred_acc_test_tmp$summary)
		pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp$summary)	
		pred_acc_test_detail = rbind(pred_acc_test_detail, pred_acc_test_tmp$prec_acc$partial_R2)
	}

	pred_acc_test = pred_acc_test[order(pred_acc_test$partial_R2, decreasing=T),]
	
	rownames(pred_acc_test_detail) = pgs_list
	colnames(pred_acc_test_detail) = paste0("rep", 1:50)
	pred_acc_test_detail = t(pred_acc_test_detail)
	pred_acc_test_detail = as.data.frame(pred_acc_test_detail)


	
	return(list(pred_acc_test=pred_acc_test, pred_acc_test_detail=pred_acc_test_detail))
}


# PRSmix

The *PRSmix* package can be used to benchmark PRS developed on an independent cohort and compute the linear combination of these score to obtained an improved prediction accuracy
# INSTALLATION
To use PRSmix:
```
install.packages("devtools")
library(devtools)
devtools::install_github("buutrg/PRSmix")
library(PRSmix)
```

# INTRODUCTION
We demonstrate the usage of PRSmix with PGS obtained from PGS catalog and evaluated on an independent cohort
- Preprocess PGS Catalog scores
- Calculate PRS with all PGS Catalog Scores
- Evaluate PRSs
- Linear combination: trait-specific (PRSmix) and cross-trait (PRSmix+)

# PREPROCESS PGS CATALOG
** Set up parameters and directory of files ** 
```
library(PRSmix)
library(foreach)
library(doParallel)

basic_data_file = "~/example/basic_data.csv"
metascore = "~/example/pgs_all_metadata_scores.csv"
score_pref = "~/data/plink_wgs"
phenofile = "~/data/phenotypes/CAD.csv"
pheno_name = "CAD_case" # Name of Phenotype in plain text
isbinary = T
pgslist = "cad_list.txt"
covar_list = c("age", "sex", paste0("PC", 1:10))
ncores = 5 # the number of CPUs process
IID_pheno = "person_id"
liabilityR2 = T
read_pred_training = F
read_pred_testing = F
```

** Extract PGS IDs from PGS catalog **

```
pgs_list_df = extract_PGSid(data = "~/data/pgs_all_metadata_scores.csv", trait=trait, efo=efo)
nrow(pgs_list_df)

fwrite(data.frame(pgs_list_df), paste0(trait, "_list.txt"), row.names=F, col.names=F, quote=F, sep="\t")

** Harmonize effect allele to ALT allele **
for (part in 0:5) {
	
	pgs_list = paste0("~/data/allPGSid.txt0", part)
	out = paste0("~/data/plink_wgs", part, ".txt")
	
	harmonize_snpeffect_toALT(
		freq_file = freq_file,
		pgs_folder = pgs_folder,
		pgs_list = pgs_list,
		out = out
	)
}

** Compute PRS ** 

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)
clusterEvalQ(cl , c(library(data.table),library(foreach)))
panel_harmonized_list = foreach(part=0:5, .export='fread') %dopar% {
	file = paste0(score_pref, anc, "_part", part, ".sscore")
	if (file.exists(file)) return(0)
	
	compute_PRS(
		plink_exe = "~/tools/plink2",
		geno = paste0("~/data/", score_pref, anc),
		ref_file = paste0("snp_weight_allpgs_", trait, "Ref_", anc, "_part", part, ".txt"),
		out = paste0(score_pref, anc, "_part", part)
		)
}
stopCluster(cl)

```

# Perform PRSmix and PRSmix+
```
score_files_list = list.files("~/data/")
score_files_list = score_files_list[grep(anc, score_files_list)]
score_files_list = score_files_list[which(endsWith(score_files_list, "sscore"))]
score_files_list = score_files_list[which(startsWith(score_files_list, score_pref))]
score_files_list = paste0("~/data/", score_files_list)
score_files_list


print(pgslist)

power_thres_list = c(0.95)
pval_thres_list = c(0.0000183)

out1 = paste0(out, "_liab")

combine_PGS(
	trait = trait,
	pgslist = pgslist,
	pheno_name = pheno_name,
	isbinary = isbinary,
	score_files_list = score_files_list,
	basic_data_file = basic_data_file,
	metascore = metascore,
	phenofile = phenofile ,
	score_pref = score_pref,
	ncores = ncores,
	IID_pheno = IID_pheno,
	liabilityR2 = liabilityR2,
	train_size_list = train_size_list,
	covar_list = covar_list,
	power_thres_list = power_thres_list,
	pval_thres_list = pval_thres_list,
	read_pred_training = read_pred_training,
	read_pred_testing = read_pred_testing,
	out = out1
)

```

The output of the combination framework contains several files:
- The case counts (for binary trait), 
- The dataframe of training and testing sample split from the main dataframe, 
- The prediction accuracy for each PRS in the training and testing set, 
- The prediction accuracy assessed in the testing set of the best PRS selected from the training set,
- the AUC of the NULL model of only covariates, the best PGS, PRSmix and PRSmix+ (adjusted for covariates), 
- Odds Ratio of the best PGS, PRSmix and PRSmix+ (adjusted for covariates), 
- The mixing weights of the scores used in combination, 
- The adjusted SNP effects to estimate PRSmix and PRSmix+ (if is_extract_adjSNPeff=TRUE)


# References
Truong, B., Hull, L.H., Ruan, Y., Huang, Q.Q., Hornsby, W., Martin, H., Heel, David v., Wang, Y., Martin, A.R.,, Lee, S.H., Natarajan, P. (2023) Integrative polygenic risk scores improve the prediction accuracy of complex traits and diseases

# Contact information
Please contact Buu Truong (btruong@broadinstitute.org) or Pradeep Natarajan (pradeep@broadinstitute.org) if you have any queries.
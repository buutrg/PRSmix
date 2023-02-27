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
- phenofile: Directory to the phenotype file
- basic_data_file: Directory to file with covariate information (age,sex,PC1..10)
- score_files_list: A vector contain directories of the PGS to be read
- pgslist: A vector of trait specific PGSs to combine
- pheno_name: Column name of the phenotype in phenofile
- isbinary: True if this is binary
- score_pref: Prefix of score files
- out: Prefix of output
- metascore: Meta-information from PGS Catalog contain PGS id and trait names. Must contains information for ALL the scores (DEFAULT = NULL)
- liabilityR2: TRUE if liability R2 should be reported (DEFAULT = FALSE)
- IID_pheno: Column name of IID of phenotype file (e.g IID, person_id)
- covar_list: A vector of of covariates, must exists as columns in basic_data_file (DEFAULT = age, sex, PC1..10))
- ncores: Number of CPU cores for parallel processing (DEFAULT = 1)
- is_extract_adjSNPeff: TRUE if extract adjSNPeff, FALSE if only calculate the combined PRS. May consume extended memory (DEFAULT = FALSE)
- snp_eff_files_list: The vector of SNP effect sizes used to compute original PRSs (DEFAULT = FALSE)
- train_size_list: A vector of training sample sizes. If NULL, all 80% of the samples will be used (DEFAULT = NULL)
- power_thres_list: A vector of power thresholds to select scores (DEFAULT = 0.95)
- pval_thres_list: A vector of P-value thresholds to select scores (DEFAULT = 0.05)
- read_pred_training: TRUE if PRSs were assessed in the training set was already run and can be read from file (DEFAULT = FALSE)
- read_pred_testing: TRUE if PRSs were assessed in the testing set was already run and can be read from file (DEFAULT = FALSE)
```

```
combine_PGS(
	phenofile,
	basic_data_file,
	score_files_list,
	pgslist,
	pheno_name,
	isbinary,
	score_pref,
	out,
	metascore,
	liabilityR2,
	IID_pheno,
	covar_list,
	ncores,
	is_extract_adjSNPeff,
	snp_eff_files_list,
	train_size_list,
	power_thres_list,
	pval_thres_list,
	read_pred_training,
	read_pred_testing
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
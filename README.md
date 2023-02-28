# PRSmix

The *PRSmix* package can be used to benchmark PRSs and compute the linear combination of trait-specific scores (PRSmix) and cross-trait scores (PRSmix+) to improve prediction accuracy.

# INSTALLATION
To use PRSmix:
```
install.packages("devtools")
library(devtools)
devtools::install_github("buutrg/PRSmix")
library(PRSmix)
```

# MANUAL
We demonstrate the usage of PRSmix with PGS obtained from PGS catalog and evaluated on an independent cohort
- Harmonize SNP effect sizes corresponding to Alternative allele in the target cohort
- Calculate PRS with all scores
- Evaluate PRSs and performed linear combination: trait-specific (PRSmix) and cross-trait (PRSmix+)

## Harmonize per-allele effect sizes to the effects of alternative allele in the target cohort

The *harmonize_snpeffect_toALT* function:
```
- ref_file: Reference file contain SNP ID (ID), reference allele (REF) and alternative allele (ALT) columns (e.g allele frequency output --freq from PLINK2)
- pgs_folder: Directory to folder contain each PGS per-allele SNP effect sizes ending with .txt
- pgs_list: File contain suffixes of file names (don't include suffix .txt) of single PGS on each line. The files must exist in the pgs_folder folder
- out: Filename of the output for the weight file
```

For example:

The *ref_file* contains 3 columns: **ID**, **REF** and **ALT** (This could be called via PLINK2 to calculate the allele frequencies with `--freq`):

| ID | ALT | REF |
| --- | --- | --- |
| rs1 | A | G |
| rs2 | T | G |
| rs3 | G | A |

In the *pgs_folder* directory (i.e. *~/example/allPGScatalog/*), there are files of per-allele SNP effect sizes: *PGS000001.txt*, *PGS000002.txt*. Each of the PGS will contains 3 columns: **SNP**, **A1**, **BETA** represent **SNP ID**, **Effect allele** and **Effect size**, respectively.

I.e. The SNP effect size file *PGS000001.txt* contains:

| SNP | A1 | BETA |
| --- | --- | --- |
| rs1 | A | 0.01 |
| rs2 | T | 0.02 |
| rs3 | G | 0.03 |

The SNP effect size file *PGS000002.txt* contains:

| SNP | A1 | BETA |
| --- | --- | --- |
| rs1 | A | 0 |
| rs2 | T | 0.03 |

The *pgs_list* (i.e. *~/example/allscoresID.txt*) file will contain (please don't include suffix .txt as this will be the name of the computed PRSs below and our package can automatically add this suffix):

```
PGS000001
PGS000002
```

Then, to harmonize SNP effect sizes:
```
harmonize_snpeffect_toALT(
	ref_file = "~/example/geno.afreq", 
	pgs_folder = "~/example/allPGScatalog/",
	pgs_list = "~/example/allscoresID.txt",
	out = "~/example/weights.txt"
)
```

The output file will contains SNP ID, A1, A2, and columns of SNP effect sizes harmonized to the A1 (alternative) allele. For example:


| SNP | A1 | A2 | PGS000001 | PGS000002 | ... |
| --- | --- | --- | --- | --- | --- |
| rs1 | A | G | 0.01 | 0 | ... |
| rs2 | T | G | 0.02 | 0.03 | ... |
| rs3 | G | A | 0.03 | 0 | ... |


## Compute PRSs for all scores

The *compute_PRS* function: 
```
- geno: Prefix of genotype file in plink format (bed/bim/fam).
- weight_file: The per-allele SNP effect output from harmonize_snpeffect_toALT function above.
- start_col: Index of the starting column of SNP effect sizes to estimate PRS in the weight file (DEFAULT = 4).
- out: Name of output file, suffix *sscore* from PLINK2 will be added.
```

Then, to compute PRSs:
```
compute_PRS(
	geno = "~/example/geno",
	weight_file = "~/example/weights.txt",
	starts_col = 4,
	out = "~/example/pgs"
)
```

## Perform linear combination: trait-specific (PRSmix) and cross-trait (PRSmix+)

The *combine_PGS* function:
```
- pheno_file: Directory to the phenotype file
- covariate_file: Directory to file with covariate information (age,sex,PC1..10)
- score_files_list: A vector contains directories of the PGSs to be read
- trait_specific_score_file: A filename contain PGS IDs of trait-specific to combine (PRSmix), one score per line
- pheno_name: Column name of the phenotype in phenofile
- isbinary: TRUE if this is binary
- out: Prefix of output
- liabilityR2: TRUE if liability R2 should be reported, otherwise partial R2 (for continuous traits) or Nagelkerke R2 (for binary traits) will be reported (DEFAULT = FALSE)
- IID_pheno: Column name of IID of phenotype file (e.g IID, person_id)
- covar_list: A vector of of covariates, must exists as columns in covariate_file (DEFAULT = age, sex, PC1..10))
- ncores: Number of CPU cores for parallel processing (DEFAULT = 1)
- is_extract_adjSNPeff: TRUE if extract adjusted SNP effects from PRSmix and PRSmix+, FALSE if only calculate the combined PRS as linear combination of PRS x mixing weights. May consume extended memory (DEFAULT = FALSE)
- original_beta_files_list: The vector contains directories to SNP effect sizes used to compute original PRSs (as weight_file argument from compute PRS above) (DEFAULT = NULL)
- train_size_list: A vector of training sample sizes. If NULL, a random 80% of the samples will be used (DEFAULT = NULL)
- power_thres_list: A vector of power thresholds to select scores (DEFAULT = 0.95)
- pval_thres_list: A vector of P-value thresholds to select scores (DEFAULT = 0.05)
- read_pred_training: TRUE if PRSs were assessed in the training set was already run and can be read from file (DEFAULT = FALSE)
- read_pred_testing: TRUE if PRSs were assessed in the testing set was already run and can be read from file (DEFAULT = FALSE)
```

The output of the combination framework contains several files:
- The number of cases and controls (for binary trait). Filename ends with `_case_counts.txt`
- The dataframe of training and testing sample split from the main dataframe. Filename ends with `_train_df.txt` and `_test_df.txt`
- The prediction accuracy in the training set for all PRSs. Filename ends with `_train_allPRS.txt`
- The prediction accuracy in the testing set for trait-specific PRSs. Filename ends with `_test_allPRS.txt`
- The prediction accuracy assessed in the testing set of the best PRS selected from the training set. Filename ends with `_best_acc.txt`
- The AUC (for binary trait) of the NULL model of only covariates, the best PGS, PRSmix and PRSmix+ (adjusted for covariates). Filename ends with `_auc_NULL.txt`, `_auc_BestPGS.txt`, `_auc_PRSmix.txt` and `_auc_PRSmixPlus.txt`
- Odds Ratio of the best PGS, PRSmix and PRSmix+ (for binary trait) (adjusted for covariates). Filename ends with  `_OR_bestPGS.txt`, `_OR_PRSmix.txt` and `_OR_PRSmixPlus.txt`
- The mixing weights of the scores used in combination. Filename ends with `_weight_PRSmix.txt` and `_weight_PRSmixPlus.txt`
- The adjusted SNP effects to estimate PRSmix and PRSmix+ (if is_extract_adjSNPeff=TRUE). Filename ends with `_adjSNPeff_PRSmix.txt` and `_adjSNPeff_PRSmixPlus.txt`
- The prediction accuracy in partial R2 or liability R2 or Nagelkerke R2 of PRSmix and PRSmix+. Filename ends with `_test_summary_traitPRS_withPRSmix.txt` and `_test_summary_traitPRS_withPRSmixPlus.txt`


For example, 

The *pheno_file* would be formatted as:

| FID | IID | CAD |
| --- | --- | --- |
| 1 | 1 | 1 |
| 2 | 2 | 0 |


Therefore,
``` 
pheno_name = "CAD"
isbinary = TRUE
liabilityR2 = TRUE
IID_pheno = "IID"
```
---
*covariate_file* as:

| FID | IID | sex | age | PC1 | ... | PC10 |
| --- | --- | --- | --- | --- | --- | --- | 
| 1 | 1 | 1 | 40 | 0.02 | ... | 0.03 |
| 2 | 2 | 0 | 50 | 0.01 | ... | 0.04 |

Therefore, 
```
covar_list = c("sex", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
```
---
*trait_specific_score_file* as (i.e. PGS IDs for *EFO_0001645* of coronary artery disease from PGS Catalog):
```
PGS000962
PGS000116
PGS000011
PGS000012
PGS000013
PGS000018
PGS000019
PGS000058
PGS000296
PGS000337
```
---
Other parameters could be:
```
ncores = 4
is_extract_adjSNPeff = TRUE
original_beta_files_list = "~/example/weights.txt"
train_size_list = NULL # randomly split 80% of the data as the training set
power_thres_list = 0.95
pval_thres_list = 0.05/2600 # P-value after Bonferroni correction
read_pred_training = FALSE
read_pred_testing = FALSE
```
---
To run the linear combination (PRSmix and PRSmix+):
```
combine_PGS(
	pheno_file = "~/example/phenotype.txt",
	covariate_file = "~/example/covariate.txt",
	score_files_list = "~/example/pgs.sscore",
	trait_specific_score_file = "cad_list.txt",
	pheno_name = "CAD",
	isbinary = TRUE,
	out = "CAD_test_",
	liabilityR2 = TRUE,
	IID_pheno = "IID",
	covar_list = c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
	ncores = 4,
	is_extract_adjSNPeff = TRUE,
	original_beta_files_list = "~/example/weights.txt",
	train_size_list = NULL,
	power_thres_list = 0.95,
	pval_thres_list = 0.05/2600,
	read_pred_training = FALSE,
	read_pred_testing = FALSE
)
```
# References
Truong, B., Hull, L. E., Ruan, Y., Huang, Q. Q., Hornsby, W., Martin, H. C., van Heel, D. A., Wang, Y., Martin, A. R., Lee, H. and Natarajan, P. 2023. "Integrative Polygenic Risk Score Improves The Prediction Accuracy Of Complex Traits And Diseases". *doi:10.1101/2023.02.21.23286110*.

# Contact information
Please contact Buu Truong (btruong@broadinstitute.org) or Pradeep Natarajan (pradeep@broadinstitute.org) if you have any queries.
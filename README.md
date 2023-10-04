# PRSmix

The *PRSmix* package can be used to benchmark PRSs and compute the linear combination of trait-specific scores (PRSmix) and cross-trait scores (PRSmix+) to improve prediction accuracy.

# INSTALLATION
To use PRSmix:
```
install.packages("devtools")
devtools::install_github("buutrg/PRSmix")
library(PRSmix)
```

# MANUAL
We demonstrate the usage of PRSmix with PGS obtained from (but not limited to) PGS catalog and evaluated on an independent cohort. The steps are:
1) Harmonize SNP effect sizes corresponding to Alternative allele in the target cohort (with `harmonize_snpeffect_toALT` function)
2) Calculate PRS with all scores (with `compute_PRS` function)
3) Evaluate PRSs and performed linear combination: trait-specific (PRSmix) and cross-trait (PRSmix+) (with `combine_PRS` function)

**BONUS**:
- If you want to evaluate a single score, you can use the `eval_single_PRS` function (described below).

**NOTE**: 
- If you already have the PRSs estimated in the target cohort (e.g. similar to [plink2 format](https://www.cog-genomics.org/plink/2.0/score)) with [`_SUM` columns](https://www.cog-genomics.org/plink/2.0/formats#sscore) (via `--score <your SNP effect file> cols=+scoresums no-mean-imputation`) and want to benchmark and combine scores, you can directly go to step 3 (evaluate and perform linear combination of the scores) with the `combine_PRS` function. If you want to return the adjusted SNP effect sizes, you will need the harmonized SNP effect sizes to the alternative allele of the target cohort (Step 1).

## Harmonize per-allele effect sizes to the effects of alternative allele in the target cohort

The *harmonize_snpeffect_toALT* function:

| Argument | Default | Description |
| --- | --- | --- |
| `ref_file` | | Reference file contain SNP ID (ID), reference allele (REF) and alternative allele (ALT) columns (e.g allele frequency output --freq from PLINK2) |
| `pgs_folder` | | Directory to folder contain each PGS per-allele SNP effect sizes ending with .txt |
| `pgs_list` | | File contain suffixes of file names (don't include suffix .txt) of single PGS on each line. The files must exist in the pgs_folder folder |
| `out` | | Filename of the output for the weight file |
| `isheader` | TRUE | TRUE if the weight files contain header, otherwise the file would be read as SNP IDs, effect allele and effect size for column 1, 2, 3, respectively |
| `snp_col` | SNP | Column name of SNP IDs. Can be a number of column index with isheader=FALSE |
| `a1_col` | A1 | Column name of effect allele. Can be a number of column index with isheader=FALSE |
| `beta_col` | BETA | Column name of effect size. Can be a number of column index with isheader=FALSE |
| `ncores` | 1 | Number of cores for parallel processing |
| `chunk_size` | 1 | Number of scores to process each chunk |

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

Optional: To split all scores into smaller chunks:
E.g. for 20 scores/chunk and 16 cores at once:
```
chunk_size = 20
ncores = 16
```

Then, to harmonize SNP effect sizes:
```
harmonize_snpeffect_toALT(
	ref_file = "~/example/geno.afreq", 
	pgs_folder = "~/example/allPGScatalog/",
	pgs_list = "~/example/allscoresID.txt",
	snp_col = 1,
	a1_col = 2,
	beta_col = 3,
	isheader = F,
	chunk_size = 20,
	ncores = 16,
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

| Argument | Default | Description |
| --- | --- | --- |
| `geno` | | Prefix of genotype file in plink format (bed/bim/fam) |
| `weight_file` | | The per-allele SNP effect output from harmonize_snpeffect_toALT function above |
| `start_col` | 4 | Index of the starting column of SNP effect sizes to estimate PRS in the weight file |
| `out` | | Name of output file, suffix *sscore* from PLINK2 will be added |


Then, to compute PRSs:
```
compute_PRS(
	geno = "~/example/geno",
	weight_file = "~/example/weights.txt",
	start_col = 4,
	out = "~/example/pgs"
)
```

The output contains multiple column names ending with `_SUM` (and `_AVG`) represent for PGS for each set of SNP effect.

| FID | IID | ALLELE_CT | NAMED_ALLELE_DOSAGE_SUM | PGS000001_AVG | PGS000001_SUM | PGS000002_AVG | PGS000002_SUM | ... |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |

See more at: https://www.cog-genomics.org/plink/2.0/formats#sscore

## Perform linear combination: trait-specific (PRSmix) and cross-trait (PRSmix+)

The *combine_PRS* function:

| Argument | Default | Description |
| --- | --- | --- |
| `pheno_file` | | Directory to the phenotype file |
| `covariate_file` | | Directory to file with covariate information (age,sex,PC1..10) |
| `score_files_list` | | A vector contains directories of the PGSs to be read as result of the `compute_PRS` function |
| `trait_specific_score_file` | | A filename contain PGS IDs of trait-specific to combine (PRSmix), one score per line |
| `pheno_name` | | Column name of the phenotype in phenofile |
| `isbinary` | | TRUE if the phenotype is a binary trait |
| `out` | | Prefix of output |
| `liabilityR2` | FALSE | TRUE if liability R2 should be reported, otherwise partial R2 (for continuous traits) or Nagelkerke R2 (for binary traits) will be reported |
| `IID_pheno` | IID | Column name of IID of phenotype file (e.g IID, person_id) |
| `covar_list` | `c("age", "sex", paste0("PC", 1:10))` | A vector of of covariates, must exists as columns in covariate_file |
| `ncores` | 1 | Number of CPU cores for parallel processing |
| `is_extract_adjSNPeff` | FALSE | TRUE if extract adjusted SNP effects from PRSmix and PRSmix+, FALSE if only calculate the combined PRS as linear combination of PRS x mixing weights. May consume extended memory |
| `original_beta_files_list` | NULL | The vector contains directories to SNP effect sizes used to compute original PRSs (as weight_file argument from compute PRS above) |
| `train_size_list` | NULL | A vector of training sample sizes. If NULL, a random 80% of the samples will be used |
| `power_thres_list` | 0.95 | A vector of power thresholds to select scores |
| `pval_thres_list` | 0.05 | A vector of P-value thresholds to select scores |
| `read_pred_training` | FALSE | TRUE if PRSs were assessed in the training set was already run and can be read from file |
| `read_pred_testing` | FALSE | TRUE if PRSs were assessed in the testing set was already run and can be read from file |


The output of the combination framework contains several files:

| Description | Suffix |
| --- | --- |
| The number of cases and controls (for binary trait) | `_case_counts.txt` |
| The dataframe of training and testing sample split from the main dataframe | `_train_df.txt` and `_test_df.txt` |
| The prediction accuracy in the training set for all PRSs | `_train_allPRS.txt` |
| The prediction accuracy in the testing set for trait-specific PRSs | `_test_allPRS.txt` |
| The prediction accuracy assessed in the testing set of the best PRS selected from the training set | `_best_acc.txt` |
| The AUC (for binary trait) of the NULL model of only covariates, the best PGS, PRSmix and PRSmix+ (adjusted for covariates) | `_auc_NULL.txt`, `_auc_BestPGS.txt`, `_auc_PRSmix.txt` and `_auc_PRSmixPlus.txt` |
| Odds Ratio of the best PGS, PRSmix and PRSmix+ (for binary trait) (adjusted for covariates) |  `_OR_bestPGS.txt`, `_OR_PRSmix.txt` and `_OR_PRSmixPlus.txt` |
| The mixing weights of the scores used in combination | `_weight_PRSmix.txt` and `_weight_PRSmixPlus.txt` |
| The adjusted SNP effects to estimate PRSmix and PRSmix+ (if is_extract_adjSNPeff=TRUE) | `_adjSNPeff_PRSmix.txt` and `_adjSNPeff_PRSmixPlus.txt` |
| The prediction accuracy of PRSmix and PRSmix+ | `_test_summary_traitPRS_withPRSmix.txt` and `_test_summary_traitPRS_withPRSmixPlus.txt` |


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
liabilityR2 = TRUE # FALSE if Nagelkerke R2 should be reported
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
combine_PRS(
	pheno_file = "~/example/phenotype.txt",
	covariate_file = "~/example/covariate.txt",
	score_files_list = c("~/example/pgs.sscore"),
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

# BONUS

The *eval_single_PRS* function can be used to evaluate a single score and return prediction accuracy R2 as partial R2, liability R2 or Nagelkerke R2:

| Argument | Default | Description |
| --- | --- | --- |
| `data_df` | | Data to assess prediction accuracy which at least contains the phenotype, covariates and a PRS |
| `pheno` | | Name of phenotype column |
| `prs_name` | | PGS list of the trait, must exist in the column name of data_df |
| `covar_list` | | Array of covariates, must exist in the column name of data_df |
| `isbinary` | FALSE | TRUE if the phenotype is a binary trait |
| `liabilityR2` | FALSE | TRUE if liability R2 should be reported, otherwise partial R2 (for continuous traits) or Nagelkerke R2 (for binary traits) will be reported |

For example: `data_df` can be formatted as:

| CAD | PRS_example | sex | age | PC1 | ... | PC10 |
| --- | --- | --- | --- | --- | --- | --- | 
| 1 | 0.02 | 1 | 40 | 0.02 | ... | 0.03 |
| 0 | 0.03 | 0 | 50 | 0.01 | ... | 0.04 |

Therefore, to evaluate `PRS_example`:
```
eval_single_PRS(
	data_df,
	pheno = "CAD",
	prs_name = "PRS_example",
	covar_list = c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
	isbinary = TRUE,
	liabilityR2 = TRUE
)
```

The output of the function will contain the R2, standard error, 95% lower and upper bound, P-value and power

# References
Truong, B., Hull, L. E., Ruan, Y., Huang, Q. Q., Hornsby, W., Martin, H. C., van Heel, D. A., Wang, Y., Martin, A. R., Lee, H. and Natarajan, P. 2023. "Integrative Polygenic Risk Score Improves The Prediction Accuracy Of Complex Traits And Diseases". *doi:10.1101/2023.02.21.23286110*.

# Contact information
Please contact Buu Truong (btruong@broadinstitute.org) or Pradeep Natarajan (pradeep@broadinstitute.org) if you have any queries.

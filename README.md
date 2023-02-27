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
- Harmonize to Alternative allele in the target cohort
- Calculate PRS with all PGS Catalog Scores
- Evaluate PRSs and performed linear combination: trait-specific (PRSmix) and cross-trait (PRSmix+)

**Harmonize to Alternative allele in the target cohort**

*harmonize_snpeffect_toALT* function:
```
- ref_file: Reference file contain SNP ID (ID), reference allele (REF) and alternative allele (ALT) columns (e.g allele frequency output --freq from PLINK2)
- pgs_folder: Directory to folder contain each PGS per-allele SNP effect sizes ending with .txt
- pgs_list: File contain suffixes of file names (ends with .txt) of single PGS on each line. The files must exist in the pgs_folder folder
- out: Filename of the output for the weight file
```

For example, in the directory *~/allPGScatalog/*, there are per-allele SNP effect sizes: PGS000001.txt PGS000002.txt. The pgs_list *~/example/allscoresID.txt* file will contain:
```
PGS000001.txt
PGS000002.txt
```

Each of the PGS will contains 3 columns: SNP, A1, BETA represent SNP ID, Effect allele and effect size

```
harmonize_snpeffect_toALT(
	ref_file = "~/example/geno.freq", 
	pgs_folder = "~/allPGScatalog/",
	pgs_list = "~/example/allscoresID.txt",
	out = "~/example/weights.txt"
)

```
**Compute PRS with all PGS Catalog Scores**

*compute_PRS* function: 
```
- geno: Prefix of genotype file in plink format (bed/bim/fam).
- weight_file: The per-allele SNP effect output from *harmonize_snpeffect_toALT*
- start_col: Index of starting column to estimate PRS in the reference file (DEFAULT = 4).
- out: Name of output file, suffix *sscore* from PLINK2 will be added.
```

```
compute_PRS(
	geno = "~/example/geno",
	weight_file = "~/example/weights.txt",
	starts_col = 4,
	out = "~/example/pgs"
)
```

**Perform linear combination: PRSmix and PRSmix+**

```
- phenofile: Directory to the phenotype file
- basic_data_file: Directory to file with covariate information (age,sex,PC1..10)
- score_files_list: A vector contain directories of the PGS to be read
- trait_specific_list: A vector of trait-specific PGS to combine
- pheno_name: Column name of the phenotype in phenofile
- isbinary: True if this is binary
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
	phenofile = "~/example/phenotype.txt",
	basic_data_file = "~/example/covariate.txt",
	score_files_list = ""~/example/pgs.sscore",
	trait_specific_list = "cad_list.txt",
	pheno_name = "CAD",
	isbinary = TRUE,
	out = "CAD_test_",
	metascore,
	liabilityR2 = TRUE,
	IID_pheno = "IID",
	covar_list = c("age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
	ncores = 4,
	is_extract_adjSNPeff = TRUE,
	snp_eff_files_list = "~/example/weights.txt",
	train_size_list = NULL,
	power_thres_list = 0.95,
	pval_thres_list = 0.05/2600,
	read_pred_training = FALSE,
	read_pred_testing = FALSE
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
Truong, B., Hull, L.H., Ruan, Y., Huang, Q.Q., Hornsby, W., Martin, H., Heel, David v., Wang, Y., Martin, A.R.,, Lee, S.H., Natarajan, P. (2023) Integrative polygenic risk scores improve the prediction accuracy of complex traits and diseases. *doi: https://doi.org/10.1101/2023.02.21.23286110*


# Contact information
Please contact Buu Truong (btruong@broadinstitute.org) or Pradeep Natarajan (pradeep@broadinstitute.org) if you have any queries.
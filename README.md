# Individual Statistics Pipeline

This repository provides a unified interface for performing SNP-level, window-based, and gene-centric association inference using knockoff-based or SCANG-based methods. It supports both **uncorrelated** and **correlated** samples, integrates SAIGE for mixed-model GWAS, and provides automated FDR control.

## Supported Methods

| Method              | Input Data            | Sample Type  | Description                                       |
| ------------------- | --------------------- | ------------ | ------------------------------------------------- |
| **KnockoffScreen**  | SNP genotypes         | Uncorrelated | SNP-level and sliding-window inference            |
| **GeneScan3DKnock** | SNP genotypes         | Uncorrelated | Gene-centric inference via multiscale aggregation |
| **BIGKnock**        | SNP genotypes (&GRM)  | Correlated   | Gene-centric inference with GLMM models           |

---

# Installation

## 1. Install SAIGE (with conda environment)

```bash
conda env create -f environment-RSAIGE.yml
conda activate RSAIGE
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
```

## 2. Install Dependencies

```bash
Rscript install_packages.R
```

---

# Input Requirements

### Required Inputs

| Argument                | Description                                                     |
| ----------------------- | --------------------------------------------------------------- |
| `--pheno_file`          | Phenotype file (CSV/TSV) containing ID column                   |
| `--geno_file`           | Genotype file prefix in PLINK bed/bim/fam format                |
| `--phenotype`           | Column name of phenotype                                        |
| `--genome_build`        | "hg19" or "hg38"                                                |
| `--sample_uncorrelated` | TRUE/FALSE                                                      |

### Optional Inputs

| Argument                    | Description                                                                                                                                | Default             |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ------------------- |
| `--sliding_window_length`   | Comma-separated list of window sizes (in bp) for window-based inference. **Only used when** `--test_type = "Single_Window"`. | `"1000,5000,10000"` |
| `--M`                       | Number of knockoff copies.                | `5`                 |
| `--geno_missing_imputation` | Missing genotype imputation method (`fixed` or `mean`).                                                                                    | `"fixed"`           |
| `--plink_path`              | Path to PLINK executable.                                                                                                              | `"plink"`           |
| `--genome_build`            | Genome build (`hg19`, `hg38`) used for annotation and window mapping.                                                                      | `"hg19"`            |
| `--sample_uncorrelated`     | Whether samples are uncorrelated (TRUE/FALSE).               | `TRUE`              |
| `--fdr`                     | Target false discovery rate.                                                                                                               | `0.1`               |
| `--grm_file`                | GRM matrix (`.grm` + `.grm.id`) **required only for correlated sample methods**, e.g., **BIGKnock**. Ignored otherwise.                    | `NULL`              |
| `--pheno_id`                | Column name of sample ID in phenotype file. Required when phenotype table contains ID-like columns.           | `NULL`              |
| `--covariates`              | Comma-separated covariate names. Optional for all methods; used only when covariates are included in association models.                   | `NULL`              |
| `--user_cores`              | Number of CPU threads used.                                                                                                                | `1`                 |


---

# Pipeline Usage

## SNP-level and window-based inference (Uncorrelated Samples)

```bash
Rscript pipeline.R \
  --outdir result/ \
  --test_type Single_Window \
  --pheno_file $pheno_file \
  --grm_file $grm_file \
  --geno_file $geno_file \
  --phenotype Y \
  --pheno_id id \
  --covariates X1 \
  --user_cores 4 \
  --sliding_window_length 1000,5000,10000 \
  --geno_missing_imputation fixed \
  --plink_path plink \
  --M 5 \
  --genome_build hg19 \
  --sample_uncorrelated TRUE \
  --fdr 0.1
```

---

## Gene-Centric Inference (Uncorrelated Samples)

```bash
Rscript pipeline.R \
  --outdir result_gene/ \
  --test_type Gene_Centric \
  --pheno_file $pheno_file \
  --geno_file $geno_file \
  --phenotype Y \
  --pheno_id id \
  --covariates X1,X2 \
  --sliding_window_length 1000,5000,10000 \
  --genome_build hg19 \
  --sample_uncorrelated TRUE \
  --M 5 \
  --fdr 0.1
```

---

## Gene-Centric Inference (Correlated Samples)

If `--grm_file` is not provided, the pipeline will generate GRM using SAIGE automatically.

```bash
Rscript pipeline.R \
  --outdir result_bigknock/ \
  --test_type Gene_Centric \
  --pheno_file $pheno_file \
  --grm_file $grm_file \
  --grm_id_file $grm_id_file \
  --geno_file $geno_file \
  --phenotype Y \
  --pheno_id id \
  --covariates age,sex,PC1,PC2 \
  --genome_build hg38 \
  --sample_uncorrelated FALSE \
  --M 5 \
  --fdr 0.1
```

---

# Output Files

Each pipeline run outputs:

### For SNP-level and Window-level inference

* `Single_Window_results_chr*.csv`
* Manhattan/Q–Q plots 

### For gene-centric inference

* `Gene_results_chr*.csv`
* Manhattan/Q–Q plots 

---

# Example Directory Structure

```
project/
├── data/
│   ├── genotype.bed
│   ├── genotype.bim
│   ├── genotype.fam
│   ├── phenotype.txt
│   └── grm/
├── R/
│   ├── pipeline.R
│   └── install_packages.R
└── result/
    ├── Single_results_chr1.csv
    ├── Window_results_chr1.csv
    └── Gene_results_chr1.csv
```

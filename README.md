# Knockoff Pipeline for Individual Statistics

## Overview

This R package provides a unified pipeline for genome-wide statistical inference using knockoff-based methods.

The pipeline supports:

- SNP-level window inference (**Single_Window**)
- Gene-centric inference (**Gene_Centric**)
- Correlated and uncorrelated sample analysis
- Mixed-model GWAS integration
- Automatic FDR control
- Batch processing with checkpoint recovery

---

## Supported Methods

| Method              | Input Data            | Sample Type  | Description                                       |
| ------------------- | --------------------- | ------------ | ------------------------------------------------- |
| **KnockoffScreen**  | SNP genotypes         | Uncorrelated | SNP-level and sliding-window inference            |
| **GeneScan3DKnock** | SNP genotypes         | Uncorrelated | Gene-centric inference via multiscale aggregation |
| **BIGKnock**        | SNP genotypes (&GRM)  | Correlated   | Gene-centric inference with GLMM models           |

---

# Installation

## 1. Install System Dependencies

### Required external software

- PLINK2

Ensure PLINK path is included in system `PATH`

## 2. Install SAIGE (with conda environment)

```bash
conda env create -f environment-RSAIGE.yml
conda activate RSAIGE
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
```

## 3. Install Other Dependencies

```bash
Rscript install_packages.R
```
## 4. Install Package

```R
devtools::install_github("tianyingw/knockoff-pipeline")
```
---

# Input Requirements

### Required Inputs

| Argument                | Description                                                     |
| ----------------------- | --------------------------------------------------------------- |
| `outdir`          | Output directory                   |
| `test_type`          | "Single_Window" or "Gene_Centric"                   |
| `pheno_file`          | Phenotype file (CSV/TSV) containing ID column                   |
| `geno_file`           | Genotype file prefix in PLINK bed/bim/fam format                |
| `phenotype`           | Column name of phenotype                                        |

### Optional Inputs

| Argument                    | Description                                                                                                                                | Default             |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ | ------------------- |
| `sliding_window_length`   | Comma-separated list of window sizes (in bp) for window-based inference. **Only used when** `--test_type = "Single_Window"`. | `"1000,5000,10000"` |
| `M`                       | Number of knockoff copies.                | `5`                 |
| `geno_missing_imputation` | Missing genotype imputation method (`fixed` or `mean`).                                                                                    | `"fixed"`           |
| `plink_path`              | Path to PLINK executable.                                                                                                              | `"plink"`           |
| `genome_build`            | Genome build (`hg19`, `hg38`) used for annotation and window mapping.                                                                      | `"hg19"`            |
| `sample_uncorrelated`     | Whether samples are uncorrelated (TRUE/FALSE).               | `TRUE`              |
| `fdr`                     | Target false discovery rate.                                                                                                               | `0.1`               |
| `grm_file`                | GRM matrix (`.grm` + `.grm.id`) **required only for correlated sample methods**, e.g., **BIGKnock**. Ignored otherwise.                    | `NULL`              |
| `pheno_id`                | Column name of sample ID in phenotype file. Required when phenotype table contains ID-like columns.           | `NULL`              |
| `covariates`              | Comma-separated covariate names. Optional for all methods; used only when covariates are included in association models.                   | `NULL`              |
| `user_cores`              | Number of CPU threads used.                                                                                                                | `1`                 |


---

# Pipeline Usage

## Load package

```R
library(KnockoffPipeline)
```

## SNP-level and window-based inference (Uncorrelated Samples)

We provide a demo at `inst/examples/SNP_Window.R`

---

## Gene-Centric Inference (Uncorrelated Samples)

We provide a demo at `inst/examples/Gene_unrelated.R`

---

## Gene-Centric Inference (Correlated Samples)

If `--grm_file` is not provided, the pipeline will generate GRM using SAIGE automatically.

We provide a demo at `inst/examples/Gene_related.R`

## Checkpoint Recovery

The pipeline supports restart from intermediate results.

Set:
```R
read_mid_exist = TRUE
```
Then pipeline will detect existing mid-results and skip finished chromosomes.

## Parallel Computing

Parallelization is supported via:
```R
user_cores = N
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



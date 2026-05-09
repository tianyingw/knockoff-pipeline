# KnockoffPipeline

## Overview

KnockoffPipeline is a unified R framework for genome-wide association analysis with rigorous false discovery rate (FDR) control via knockoff-based methods.

The pipeline supports:

- SNP-level and sliding-window inference (**Single_Window**)
- Gene-centric inference with 3D enhancer information (**Gene_Centric**)
- Unrelated and related sample analysis
- Multiple phenotypes in a single run (shared knockoffs across phenotypes)
- **Knockoff persistence**: save generated knockoffs and reload them across sessions
- **Two-stage workflow**: decouple knockoff generation from association testing

---

## Supported Methods

| Method              | Input              | Sample type  | Description                                          |
|---------------------|--------------------|--------------|------------------------------------------------------|
| **KnockoffScreen**  | SNP genotypes      | Unrelated    | SNP-level and sliding-window inference               |
| **GeneScan3DKnock** | SNP genotypes      | Unrelated    | Gene-centric inference via 3D enhancer mapping       |
| **BIGKnock**        | SNP genotypes + GRM | Related     | Gene-centric inference with GLMM null model          |

---

## Installation

### conda

#### 1. create a conda environment

```bash
conda env create -f inst/conda_env/environment.yml
conda activate pipeline
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
```

#### 2. Check and install R dependencies

```bash
Rscript inst/conda_env/install_packages.R
```

#### 3. Install SAIGE

```bash
R CMD INSTALL SAIGE
```

#### 4. Install KnockoffPipeline

```R
devtools::install_github("tianyingw/knockoff-pipeline")
```

---

### docker

This image targets Linux containers.

Build the image from the package directory:

```bash
cd knockoff-pipeline
docker build -f inst/docker/Dockerfile -t knockoff-pipeline:latest .
```

If you want to force regeneration of `SAIGE/plink2_includes.a` inside the
container instead of reusing the checked-in archive:

```bash
docker build -f inst/docker/Dockerfile \
  --build-arg REBUILD_PLINK2_INCLUDES=1 \
  -t knockoff-pipeline:latest .
```

Open a shell in the container:

```bash
docker run --rm -it knockoff-pipeline:latest
```

Inside the container, the working directory is already the package root, so the
example scripts can be run directly:

```bash
Rscript inst/examples/SNP_Window.R
Rscript inst/examples/Gene_unrelated.R
Rscript inst/examples/Gene_related.R
```

You can also run them directly from the host:

```bash
docker run --rm knockoff-pipeline:latest Rscript inst/examples/SNP_Window.R
docker run --rm knockoff-pipeline:latest Rscript inst/examples/Gene_unrelated.R
docker run --rm knockoff-pipeline:latest Rscript inst/examples/Gene_related.R
```

If you want the example output on the host:

```bash
mkdir -p results
docker run --rm \
  -e KNOCKOFF_OUTDIR=/results \
  -v "$PWD/results:/results" \
  knockoff-pipeline:latest \
  Rscript inst/examples/SNP_Window.R
```

Run the pipeline on your own data:

```bash
docker run --rm \
  -v /absolute/path/to/data:/data \
  -v /absolute/path/to/results:/results \
  knockoff-pipeline:latest \
  Rscript -e "library(KnockoffPipeline); run_pipeline(outdir='/results', test_type='Single_Window', pheno_file='/data/phenotype.csv', geno_file='/data/demo', phenotype='Y', pheno_id='IID', covar_cols='X1', plink_path='plink2', genome_build='hg19', sample_uncorrelated=TRUE)"
```

Notes:

- The image installs the local `SAIGE` source tree first, then installs `KnockoffPipeline`.
- `plink2` is available inside the container from the conda environment, so use `plink_path = "plink2"` unless you intentionally mount another binary.
- The bundled example scripts accept `KNOCKOFF_OUTDIR`, `PLINK_BIN`, and `KNOCKOFF_CORES` environment variables.
- The current Dockerfile reuses the checked-in `SAIGE/plink2_includes.a` by default. Set `REBUILD_PLINK2_INCLUDES=1` if you need to rebuild it in-container.
- Build from the `knockoff-pipeline/` subdirectory, not the repository root, to avoid sending large unrelated files into the Docker build context.

### apptainer / singularity

This is the recommended container route on many HPC servers where Docker is not available.

Build the image from the package directory:

```bash
cd knockoff-pipeline
apptainer build knockoff-pipeline.sif inst/apptainer/knockoff-pipeline.def
```

If your cluster still uses the older command name, `singularity` is equivalent:

```bash
singularity build knockoff-pipeline.sif inst/apptainer/knockoff-pipeline.def
```

Open a shell:

```bash
apptainer shell knockoff-pipeline.sif
```

Run the packaged examples directly:

```bash
apptainer exec knockoff-pipeline.sif Rscript inst/examples/SNP_Window.R
apptainer exec knockoff-pipeline.sif Rscript inst/examples/Gene_unrelated.R
apptainer exec knockoff-pipeline.sif Rscript inst/examples/Gene_related.R
```

Write example output to a host directory:

```bash
mkdir -p results
apptainer exec \
  --bind "$PWD/results:/results" \
  --env KNOCKOFF_OUTDIR=/results \
  knockoff-pipeline.sif \
  Rscript inst/examples/SNP_Window.R
```

Run on your own data:

```bash
apptainer exec \
  --bind /absolute/path/to/data:/data \
  --bind /absolute/path/to/results:/results \
  knockoff-pipeline.sif \
  Rscript -e "library(KnockoffPipeline); run_pipeline(outdir='/results', test_type='Single_Window', pheno_file='/data/phenotype.csv', geno_file='/data/demo', phenotype='Y', pheno_id='IID', covar_cols='X1', plink_path='plink2', genome_build='hg19', sample_uncorrelated=TRUE)"
```

Notes:

- The definition file is [inst/apptainer/knockoff-pipeline.def](/abs/path/e:/college_material/statistics/tianying/knockoff-pipeline/knockoff-pipeline/inst/apptainer/knockoff-pipeline.def).
- It reuses the checked-in `SAIGE/plink2_includes.a`, so the current setup targets Linux `x86_64`.
- `apptainer build` often requires root, `sudo`, or a site-supported remote build workflow. On many HPC systems, users can `exec` and `shell` existing `.sif` files but cannot build them locally.
- If your cluster forbids local builds, build `knockoff-pipeline.sif` on another Linux machine first, then upload the `.sif` file to the server.


## Quick Start

```R
library(KnockoffPipeline)

run_pipeline(
  outdir     = "results/",
  test_type  = "Single_Window",   # or "Gene_Centric"
  pheno_file = "data/pheno.csv",
  geno_file  = "data/geno",       # PLINK prefix
  phenotype  = "BMI"
)
```

---

## Input Requirements

### Required

| Argument     | Description                                              |
|--------------|----------------------------------------------------------|
| `outdir`     | Output directory (created automatically if absent)       |
| `test_type`  | `"Single_Window"` or `"Gene_Centric"`                   |
| `pheno_file` | Phenotype file (CSV/TSV)                                 |
| `geno_file`  | PLINK genotype file prefix (`.bed/.bim/.fam`)            |
| `phenotype`  | Column name(s) of phenotype(s) — scalar or character vector |

### Optional

| Argument                  | Description                                                                                      | Default               |
|---------------------------|--------------------------------------------------------------------------------------------------|-----------------------|
| `pheno_id`                | Column name of sample ID in phenotype file                                                       | `NULL`                |
| `covar_cols`              | Continuous covariate column names                                                                | `NULL`                |
| `cat_covar_cols`          | Categorical covariate column names                                                               | `NULL`                |
| `sliding_window_length`   | Window size (bp) for `Single_Window` mode                                                        | `"1000,5000,10000"`   |
| `M`                       | Number of knockoff copies                                                                        | `5`                   |
| `geno_missing_imputation` | Genotype imputation method (`"fixed"` or `"mean"`)                                               | `"fixed"`             |
| `plink_path`              | Path to PLINK executable                                                                         | `"plink"`             |
| `genome_build`            | `"hg19"` or `"hg38"`                                                                             | `"hg19"`              |
| `sample_uncorrelated`     | `TRUE` = GLM null model; `FALSE` = GLMM via SAIGE                                               | `TRUE`                |
| `grm_file`                | Sparse GRM file (required only for `sample_uncorrelated = FALSE`)                                | `NULL`                |
| `fdr`                     | Target FDR level                                                                                 | `0.1`                 |
| `chromosomes`             | Autosomes to analyse                                                                             | `1:22`                |
| `user_cores`              | Number of CPU threads                                                                            | `1`                   |
| `batch_size`              | Genes per batch (Gene_Centric only)                                                              | `20`                  |
| `read_mid_exist`          | Skip chromosomes with existing intermediate files                                                | `TRUE`                |
| **`pipeline_stage`**      | `"full"`, `"stage1_knockoff"`, or `"stage2_analysis"` — see below                              | `"full"`              |
| **`save_knockoff`**       | Save generated knockoffs to disk (`TRUE`/`FALSE`/`NULL` for auto)                               | `NULL`                |
| **`knockoff_dir`**        | Directory for knockoff `.rds` files (defaults to `<outdir>/knockoffs`)                           | `NULL`                |

---

## Key Features

### Multiple Phenotypes

Pass a character vector to `phenotype`. The pipeline:

1. Removes samples missing in **any** phenotype or covariate once, producing a single consistent sample set.
2. Generates knockoffs on the first phenotype pass and **automatically saves and reloads** them for all subsequent phenotypes (since knockoffs depend only on genotype, not phenotype).
3. Writes per-phenotype results to `<outdir>/<phenotype_name>/`.

```R
run_pipeline(
  outdir    = "results/",
  test_type = "Gene_Centric",
  pheno_file = "data/pheno.csv",
  geno_file  = "data/geno",
  phenotype  = c("BMI", "LDL", "SBP"),   # three phenotypes
  pheno_id   = "IID"
)
```

### Knockoff Persistence (`save_knockoff`)

Set `save_knockoff = TRUE` to save generated knockoffs to `knockoff_dir` (one `.rds` file per LD block or per gene, under `<knockoff_dir>/chr<c>/`). On subsequent analyses with the same genotype data but a different phenotype file, set `load_knockoff = TRUE` (handled automatically when `save_knockoff = TRUE`) to skip knockoff generation entirely.

```R
# First run: generate and save
run_pipeline(..., save_knockoff = TRUE)

# Later run: load saved knockoffs
run_pipeline(..., pipeline_stage = "stage2_analysis", knockoff_dir = "results/knockoffs")
```

Each knockoff file stores the sample ID list and SNP positions. On load, the pipeline checks that:
- The **column count** (number of SNPs after QC) matches. If not, a warning is issued and knockoffs are regenerated.
- The **sample order** is aligned to the current analysis. Missing samples trigger a warning; extra saved samples are silently ignored.

### Two-Stage Workflow (`pipeline_stage`)

The pipeline can be split into two independent stages, useful when knockoff generation and association testing need to run in separate jobs (e.g., on a cluster), or when the same knockoffs will be reused with multiple phenotype files that are not yet available.

| `pipeline_stage`       | What it does                                                                 |
|------------------------|------------------------------------------------------------------------------|
| `"full"` (default)     | Complete end-to-end pipeline                                                 |
| `"stage1_knockoff"`    | Generate and save knockoffs only; write sample list; no association testing  |
| `"stage2_analysis"`    | Load saved knockoffs; fit null models; run association tests                 |

**Stage 1:**

```R
run_pipeline(
  outdir         = "results/",
  test_type      = "Gene_Centric",
  pheno_file     = "data/pheno.csv",   # used only to determine the sample set
  geno_file      = "data/geno",
  phenotype      = "BMI",
  pipeline_stage = "stage1_knockoff",
  knockoff_dir   = "results/knockoffs"
)
# Outputs:
#   results/knockoffs/chr1/gene_BRCA1_ko.rds  ...
#   results/knockoffs/knockoff_sample_list.txt
```

**Stage 2** (can use a completely different phenotype file):

```R
run_pipeline(
  outdir         = "results_AD/",
  test_type      = "Gene_Centric",
  pheno_file     = "data/pheno_AD.csv",
  geno_file      = "data/geno",
  phenotype      = "AD",
  pipeline_stage = "stage2_analysis",
  knockoff_dir   = "results/knockoffs"   # same knockoffs from stage 1
)
```

Stage 2 reads `knockoff_sample_list.txt`, computes the intersection with the current genotype/phenotype data, reindexes knockoff rows accordingly, and warns if any samples are lost in either direction.

---

## Checkpoint Recovery

The pipeline supports automatic restart from intermediate results. Set:

```R
read_mid_exist = TRUE   # (default)
```

The pipeline will detect existing per-chromosome intermediate files and skip completed chromosomes.

---

## Output Files

### Single_Window

| File                                  | Description                      |
|---------------------------------------|----------------------------------|
| `Single_Window_results.csv`           | Full results table               |
| `manhattan_plot_single.png`           | Manhattan / Q–Q plot             |
| `mid/Single_mid_results_chr*.txt`     | Per-chromosome intermediate SNP  |
| `mid/Window_mid_results_chr*.txt`     | Per-chromosome intermediate window |

### Gene_Centric

| File                                  | Description                      |
|---------------------------------------|----------------------------------|
| `GeneCentric_results.csv`             | Full results table               |
| `manhattan_plot_gene.png`             | Manhattan / Q–Q plot             |
| `mid/GeneCentric_mid_results_chr*.txt`| Per-chromosome intermediate      |

### Knockoffs (when `save_knockoff = TRUE`)

| File                                           | Description                                 |
|------------------------------------------------|---------------------------------------------|
| `knockoffs/knockoff_sample_list.txt`           | Sample IID list in knockoff row order       |
| `knockoffs/chr<c>/block_XXXX_knockoff.rds`    | Per-LD-block knockoff (Single_Window)       |
| `knockoffs/chr<c>/gene_<ID>_ko.rds`           | Per-gene knockoff, gene buffer only (Gene_Centric) |

### Multi-phenotype runs

Results for each phenotype are written to `<outdir>/<phenotype_name>/`.

---

## Examples

Demo scripts are provided under `inst/examples/`:

| Script                 | Description                                        |
|------------------------|----------------------------------------------------|
| `SNP_Window.R`         | SNP/window-level analysis, unrelated samples       |
| `Gene_unrelated.R`     | Gene-centric analysis, unrelated samples           |
| `Gene_related.R`       | Gene-centric analysis, related samples (BIGKnock)  |
<!-- | `MultiPheno.R`         | Multi-phenotype run with shared knockoffs          |
| `TwoStage.R`           | Stage 1 + Stage 2 split workflow                   | -->

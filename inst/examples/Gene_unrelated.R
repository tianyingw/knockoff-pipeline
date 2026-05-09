library(KnockoffPipeline)

example_outdir <- Sys.getenv("KNOCKOFF_OUTDIR", unset = "inst/examples/output")
plink_bin <- Sys.getenv("PLINK_BIN", unset = "plink2")
example_cores <- as.integer(Sys.getenv("KNOCKOFF_CORES", unset = "1"))

run_pipeline(
    outdir = example_outdir,
    test_type = "Gene_Centric",
    pheno_file = "inst/examples/input/phenotype.csv",
    geno_file = "inst/examples/input/demo",
    phenotype = "Y",
    pheno_id = "IID",
    covar_cols = c("X1"),
    user_cores = example_cores,
    sliding_window_length = c(1000,5000,10000),
    geno_missing_imputation = "fixed",
    plink_path = plink_bin,
    M = 5,
    genome_build = "hg19",
    sample_uncorrelated = TRUE,
    fdr = 0.1,
    chromosomes = 1
)

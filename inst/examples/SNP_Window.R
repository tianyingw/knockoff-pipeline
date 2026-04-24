library(KnockoffPipeline)

run_pipeline(
    outdir = "inst/examples/output",
    test_type = "Single_Window",
    pheno_file = "inst/examples/input/phenotype.csv",
    geno_file = "inst/examples/input/demo",
    phenotype = "Y",
    pheno_id = "IID",
    covar_cols = c("X1"),
    user_cores = 8,             
    sliding_window_length = c(1000,5000,10000),
    geno_missing_imputation = "fixed",
    plink_path = "/dssg/home/acct-mashiyang1991/mashiyang1991/yujie/ukb_data/imputed/plink2",
    M = 5,
    genome_build = "hg19",
    sample_uncorrelated = TRUE,
    fdr = 0.1,
    chromosome = 1
)
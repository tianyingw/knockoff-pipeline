library(KnockoffPipeline)

run_pipeline(
    outdir = "output",
    test_type = "Gene_Centric",
    pheno_file = "input/phenotype.csv",
    geno_file = "input/demo",
    phenotype = "Y",
    pheno_id = "IID",
    covariates = c("age","sex"),
    user_cores = 8,             
    sliding_window_length = c(1000,5000,10000),
    geno_missing_imputation = "fixed",
    plink_path = "your_path_to_plink/plink",
    M = 5,
    genome_build = "hg19",
    sample_uncorrelated = FALSE,
    fdr = 0.1,
    chromosome = 1
)
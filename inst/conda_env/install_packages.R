#!/usr/bin/env Rscript

# # install SAIGE
#install required R packages, from Finnge/SAIGE-IT
req_packages <- c("R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "roxygen2", "rversions","devtools", "SKAT", "RhpcBLASctl", "qlcMatrix", "RSQLite", "lintools")
for (pack in req_packages) {
    if(!require(pack, character.only = TRUE)) {
        install.packages(pack, repos = "https://cloud.r-project.org")
    }
}

#devtools::install_github("leeshawn/SPAtest")
devtools::install_github("leeshawn/MetaSKAT")
#devtools::install_github("leeshawn/SKAT")
devtools::install_github('chrchang/plink-ng', subdir='2.0/pgenlibr')

## Install required R packages for the pipeline (CRAN + GitHub)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
packages_to_load <- c("Matrix", "bigmemory", "CompQuadForm", "data.table","SPAtest", "irlba", 
                      "SKAT", "MASS", "WGScan", "abind",
                     "tictoc", "dplyr", "parallel", "qqman", "devtools")

missing_packages <- packages_to_load[!packages_to_load %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  message("正在安装缺失的包: ", paste(missing_packages, collapse = ", "))

  # 安装缺失的包
  for (pkg in missing_packages) {
    message("安装包: ", pkg)
    tryCatch({
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
      message("成功安装: ", pkg)
    }, error = function(e) {
      warning("安装失败: ", pkg, " - ", e$message)
    })
  }
} else {
  message("所有需要的CRAN-R包都已安装。")
}

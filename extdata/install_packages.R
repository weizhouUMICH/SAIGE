#!/usr/bin/env Rscript
#install required R packages, from Finnge/SAIGE-IT

#req_packages <- c("R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "MetaSKAT", "roxygen2", "rversions","devtools", "SKAT")
req_packages <- c("R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "roxygen2", "rversions","devtools", "SKAT", "RhpcBLASctl")
for (pack in req_packages) {
    if(!require(pack, character.only = TRUE)) {
        install.packages(pack, repos = "https://cloud.r-project.org")
    }
}

#devtools::install_github("leeshawn/SPAtest")
devtools::install_github("leeshawn/MetaSKAT")
#devtools::install_github("leeshawn/SKAT")

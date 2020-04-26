#!/usr/bin/env Rscript
packages = c("Rcpp", "RcppParallel", "parallel")
for(package in packages){
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    if (!require(package, character.only = TRUE)) stop("Load failure: ", package)
  }
}

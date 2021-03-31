#' TTIQ-CIscript.R
#' Author: Peter Ashcroft, ETH Zurich

library(parallel)
library(tidyverse)

ParLapply <- function(X, FUN, ..., PARALLEL = TRUE, SEED = NULL) {
  if (PARALLEL) {
    if (is.null(SEED)) SEED <- 111
    num.cores <- min(c(length(X), max(detectCores() - 1, 1)))
    cl <- makeCluster(num.cores, type = "FORK")
    clusterSetRNGStream(cl, SEED)
    out <- parLapply(cl, X, FUN, ...)
    stopCluster(cl)
    return(out)
  } else {
    out <- lapply(X, FUN, ...)
    return(out)
  }
}

getCI <- function(paramList, output) {
  #' For each infProf param set, we want to know which genDist params are within
  #' the joint 95% CI
  minLLH <- max(infParamsLLH$llh) + max(genParamsLLH$llh) - qchisq(0.95, df = 3+2)/2
  genParamID <- lapply(seq_len(nrow(infParamsLLH)), function(id) {
    which(genParamsLLH$llh + infParamsLLH[id,"llh"] > minLLH)
  })

  cases.outer <- ParLapply(seq_len(nrow(infParamsLLH)), function(id) {
    cases.inner <- lapply(genParamID[[id]], function(idGD) {
      #' Calculate cases for these parameter sets
      out.df <- func(id, idGD)
      return(out.df)
    })
    #' List of dataframes for each idGD: want to extract lower and upper
    df <- cases.inner[[1]]
    lower <- df[[output]]
    upper <- df[[output]]
    for (idGD in seq_along(genParamID[[id]])) {
      lower <- pmin(lower, cases.inner[[idGD]][[output]])
      upper <- pmax(upper, cases.inner[[idGD]][[output]])
    }
    #' Add to DF
    df$lower <- lower
    df$upper <- upper
    return(df[,c(names(paramList), "lower", "upper")])
  })
  #' List of dataframes for each id: want to extract lower and upper
  df <- cases.outer[[1]]
  lower <- df$lower
  upper <- df$upper
  for (id in seq_len(nrow(infParamsLLH))) {
    lower <- pmin(lower, cases.outer[[id]]$lower)
    upper <- pmax(upper, cases.outer[[id]]$upper)
  }
  #' Add to DF
  df$lower <- lower
  df$upper <- upper
  return(df)
}

checkFiles <- function(info.filename) {
  #' Check we have the information to run
  if (!file.exists(info.filename)) stop(paste("Info file", info.filename, "is missing"))
  #' Check the data directory is present
  if (!dir.exists("data/"))  stop("Data directory /data is missing")
  #' Check for the four data files which are required
  for (filename in c("data/savedDistributions.RData", "data/savedDistributionsLLH.RData",
                     "data/infParamsLLH.RData", "data/genParamsLLH.RData")) {
    if (!file.exists(filename)) stop(paste("Data file", filename, "is missing"))
  }
}

# Evaluate on Euler (HPC) ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Missing data filename")
} else if (length(args) == 1) {
  info.filename <- args[1]
  checkFiles(info.filename)
  #' Submit the jobs
  #invisible(system(paste("echo", "Submitting the jobs", sep = " ")))
  invisible(system(paste0("bsub -n 48 -R fullnode \"", paste("R --vanilla --slave < TTIQ-CIscript.R --args", info.filename, "run", sep = " "), "\"")))

} else if (length(args) == 2) {
  info.filename <- args[1]
  checkFiles(info.filename)
  #' Load the information
  load(info.filename)
  #' Load the data files
  load("data/savedDistributions.RData")
  load("data/savedDistributionsLLH.RData")
  load("data/infParamsLLH.RData")
  load("data/genParamsLLH.RData")
  #' Execute the job and save final dataframe
  out.df <- getCI(paramList = paramList, output = output)
  save(out.df, file = out.filename)
}

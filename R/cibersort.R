# Public code https://rdrr.io/github/IOBR/IOBR/src/R/CIBERSORT.R

# CIBERSORT R script v1.03 (last updated 07-10-2015)
# Note: Signature matrix construction is not currently available; use java version for full functionality.
# Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
# Requirements:
#       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#       install.packages('e1071')
#       install.pacakges('parallel')
#       install.packages('preprocessCore')
#       if preprocessCore is not available in the repositories you have selected, run the following:
#           source("http://bioconductor.org/biocLite.R")
#           biocLite("preprocessCore")
# Windows users using the R GUI may need to Run as Administrator to install or update packages.
# This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
# single-threaded in Windows.
#
# Usage:
#       Navigate to directory containing R script
#
#   In R:
#       source('CIBERSORT.R')
#       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#
#       Options:
#       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#       ii) QN = Quantile normalization of input mixture (default = TRUE)
#
# Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
# Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
# License: http://cibersort.stanford.edu/CIBERSORT_License.txt


# Core algorithm of Cibersort
#'
#' @keywords internal
#'
#' @importFrom parallel mclapply
#' @importFrom stats cor
#'
#'
CoreAlg <- function(X, y, cores = 3){


  ########################
  ## X is the data set
  ## y is labels for each row in X
  ########################



  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}

    #if(i==1){nus <- 0.997}
    #if(i==2){nus <- 0.998}
    #if(i==3){nus <- 0.999}

    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=FALSE)
    model
  }

  #Execute In a parallel way the SVM
  if(cores>1){
    if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1)
    else out <-  parallel::mclapply(1:svn_itor, res, mc.cores=cores)
  }
  else out <-  lapply(1:svn_itor, res)

  #Initiate two variables with 0
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  ##############################
  ## Here CIBERSORT starts    #
  ##############################

  t <- 1
  while(t <= svn_itor) {

    #Get the weights with a matrix multiplications between two vectors. I should get just one number (?)
    #This is done multiplying the coefficients (?) and ???

    #The support vectors
    #are the points of my dataset that lie closely to the plane that separates categories
    #The problem now is that I don't have any category (discrete variable, e.g., "sport", "cinema") but I ave continuous variable
    mySupportVectors <- out[[t]]$SV

    #My coefficients
    myCoefficients <- out[[t]]$coefs

    weights = t(myCoefficients) %*% mySupportVectors

    #set up weight/relevance on each
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)

    #This multiplies the reference profile for the correspondent weigth
    u <- sweep(X,MARGIN=2,w,'*')

    #This does the row sums
    k <- apply(u, 1, sum)

    #Don't know
    nusvm[t] <- sqrt((mean((k - y)^2))) #pitagora theorem
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  #print(mn)
  model <- out[[mn]]

  #get and normalize coefficients

  #############################################
  ## THIS IS THE SECRET OF CIBERSORT
  #############################################

  q <- t(model$coefs) %*% model$SV

  #############################################
  #############################################

  q[which(q<0)]<-0

  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}

#' @importFrom stats sd
#'
#' @keywords internal
#'
doPerm <- function(perm, X, Y, cores = 3){


  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr, cores = cores)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}

# MADE BY STEFANO TO ALLOW PARALLELISM
call_core = function(itor, Y, X, P, pval, CoreAlg){
  ##################################
  ## Analyze the first mixed sample
  ##################################

  y <- Y[,itor]

  #standardize mixture
  y <- (y - mean(y)) / sd(y)

  #run SVR core algorithm
  result <- CoreAlg(X, y, cores = 1)

  #get results
  w <- result$w
  mix_r <- result$mix_r
  mix_rmse <- result$mix_rmse

  #calculate p-value
  if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

  #print output
  c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)

}


#' @importFrom stats sd
#' @importFrom utils install.packages
#'
#' @keywords internal
#'
my_CIBERSORT <- function(Y, X, perm=0, QN=TRUE, cores = 3, exp_transform = FALSE){


  #read in data
  X <- data.matrix(X)
  Y <- data.matrix(Y)

  #order

  ###################################
  ## This is needed to make the two tables consistent in gene
  ###################################
  common_genes = intersect(rownames(X), rownames(Y))
  X <- X[common_genes,,drop=FALSE]
  Y <- Y[common_genes,,drop=FALSE]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(is.null(exp_transform)) exp_transform = max(Y) < 50
  if(exp_transform) {Y <- 2^Y}

  #quantile normalization of mixture file

  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,,drop=FALSE]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,,drop=FALSE]

  # Eliminate empty samples
  if(length(which(colSums(Y)==0))>0)
    warning(sprintf(
      "tidybulk says: the samples %s were ignored for decovolution as they have 0 counts for the deconvolution signature genes",
      colnames(Y)[colSums(Y)==0] %>% paste(collapse = ", ")
    ))
  Y=Y[,colSums(Y)>0, drop=FALSE]

  # Check if package is installed, otherwise install
  if (find.package("matrixStats", quiet = TRUE) %>% length %>% equals(0)) {
    message("tidybulk says: Installing matrixStats needed for cibersort")
    install.packages("matrixStats", repos = "https://cloud.r-project.org")
  }

  # Eliminate sd == 0
  if(length(which(matrixStats::colSds(Y)==0))>0)
    warning(sprintf(
      "tidybulk says: the samples %s were ignored for decovolution as they have standard deviation of 0 for the deconvolution signature genes",
      colnames(Y)[matrixStats::colSds(Y)==0] %>% paste(collapse = ", ")
    ))
  Y = Y[,matrixStats::colSds(Y)>0,drop=FALSE]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  # stefano write X and Y
  Y_norm <- apply(Y, 2, function(mc) (mc - mean(mc)) / sd(mc)  )

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y, cores = cores)$dist)}


  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")

  output <- matrix()
  itor <- 1
  mix <- dim(Y)[2]
  pval <- 9999

  # If not Windows
  if(Sys.info()['sysname'] == 'Windows')
  {
    while(itor <= mix){

      ##################################
      ## Analyze the first mixed sample
      ##################################


      out <- call_core(itor, Y, X, P, pval, CoreAlg)
      if(itor == 1) {output <- out}
      else {output <- rbind(output, out)}
      itor <- itor + 1

    }

  }

  # If Linux of Mac
  else {
    output <- parallel::mclapply(1:mix, call_core, Y, X, P, pval, CoreAlg, mc.cores=cores)
    output= matrix(unlist(output), nrow=length(output), byrow=TRUE)

  }

  #save results
  #write.table(rbind(header,output), file="CIBERSORT-Results.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1, drop=FALSE]
  obj <- obj[-1,, drop=FALSE]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")

  list(proportions=obj, mix = Y_norm, signatures = X)

}

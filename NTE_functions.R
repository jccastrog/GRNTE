#!/usr/bin/env Rscript
################################################################################
#Name:    NTE_functions.R
#Author:  Juan C. Castro <jcastro37@gatech.edu>
#         Diego M. Riaño P. <diriano@gmail.com>
#Update:  07-Feb-2018
#Version: 1.0.2
#License: GNU General Public License v3.0.
#===============================================================================
#
################################################################################
#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
  dir.create(personal.lib.path)

packages <- c("optparse", "entropy", "gdata")
if(any(!(packages %in% installed.packages()))){
  cat("Please wait a moment! Installing required packages ...\n")
  install.packages(packages[!(packages %in% installed.packages())],
                   quiet = T, repos="http://cran.rstudio.com/",
                   lib = personal.lib.path)
  cat("Required packages installed!\n")
}

#======= 1.0 Load packages, define functions, and initialize variables =======#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(entropy))
suppressPackageStartupMessages(library(gdata))

# 1.2 Define functions =======================================================#
#' Estimate mutual information (MI) for all pairs of variables in an expression
#' matrix
#'
#' @param matrixData A matrix with expression data where columns are genes and
#'        rows are time points
#' @param numGenes The number of genes in the matrix
#' @return mutualMat A matrix with values of mutual information for all pair of
#'        genes 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
mutualInfoEst <- function(matrixData,numGenes) {
  numBins <- ceiling(log2(nrow(matrixData)+1))
  mutualMat <- matrix(ncol=numGenes,nrow=numGenes)
  for (i in 1:numGenes) {
    for (j in 1:numGenes) {
      discretVec <- discretize2d(matrixData[,i],matrixData[,j],numBins,numBins)
      mutualMat[i,j] <- suppressWarnings(mi.empirical(discretVec,unit=c("log2")))
    }
  }
  return(mutualMat)
}
#' Estimate a null distribution of MI for variables in an expression matrix
#'
#' @param matrixData A matrix with expression data where columns are genes and
#'        rows are time points
#' @param numGenes The number of genes in the matrix
#' @param randomizations Number of times to randomize expression values
#' @return distMat A matrix with values of mutual information for a randomized
#'        expression matrix where each column is a linearized version of a random
#'        expression matrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
nullInfoDist <- function(matrixData,numGenes,randomizations){
  distMat <- matrix(ncol=(numGenes^2),nrow=randomizations)
  for (counts in 1:randomizations) {
    randomMatrixData <- matrix(ncol=ncol(matrixData),nrow=(nrow(matrixData)))
    for (i in 1:ncol(matrixData)){
      randomMatrixData[,i] <- sample(matrixData[,i])
    }
    randomMIMat <- mutualInfoEst(randomMatrixData,numGenes)
    randomMIVec <- matrix(randomMIMat,ncol=(numGenes^2))
    distMat[counts,] <- randomMIVec
  }
  return(distMat)
}
#' Calculate a score for each value of mutual information in an MI matrix
#'
#' @param initialMatrix A matrix with MI values calculated form an expression
#'        matrix
#' @param nullDistMaxtrix A matrix with null values of MI obatined with nullInfoDist
#' @param numGenes The number of genes in the matrix
#' @param randomizations Number of times to randomize expression values
#' @return pMat A matrix with pValues of for the MI values in initialMatrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
infoScore <- function(initialMatrix,nullDistMatrix,numGenes,randomizations){
  matSize <- numGenes^2
  sumDist <- randomizations
  linealIniMat <- matrix(initialMatrix,ncol=matSize)
  pVec <- c()
  for (i in 1:matSize){
    MICounts <- sum(nullDistMatrix[,i] >= linealIniMat[i])
    pScore <- MICounts/sumDist
    pVec[i] <- pScore
  }
  pMat <- matrix(pVec,ncol=numGenes,nrow=numGenes)
  return(pMat)
}
#' Calculate the MI values for an expression matrix with shifted values 
#'
#' @param matrixData A matrix with expression data where columns are genes and
#'        rows are time points
#' @param numRep The number of replicates for each time point
#' @param stepSize The size of the step to take each time.
#' @param randomizations Number of times to randomize expression values
#' @return lagMatrix An asymetric matrix with values of MI for shifted variables
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
lagMIMat <- function(matrixData,numRep, stepSize){
  numGenes <- ncol(matrixData)
  numBins <- ceiling(log2(nrow(matrixData)+1))
  step <- numRep * stepSize
  headMatrix <-matrixData[1:(nrow(matrixData)-step),]
  tailMatrix <- matrixData[(step+1):(nrow(matrixData)),]
  lagMatrix <- matrix(ncol=numGenes,nrow=numGenes)
  for (i in 1:numGenes){
    for (j in 1:numGenes){
      lagVec<- suppressWarnings(discretize2d(headMatrix[,i],tailMatrix[,j],numBins,numBins))
      lagMatrix[i,j] <- suppressWarnings(mi.empirical(lagVec,unit=c("log2")))
    }
  }
  return(lagMatrix)
}
#' Calculate the null distribution of  MI values for an expression matrix 
#' with shifted values obtained with lagMIMat
#'
#' @param matrixData A matrix with expression data where columns are genes and
#'        rows are time points
#' @param numRep The number of replicates for each time point
#' @param randomizations Number of times to randomize expression values
#' @return distMat A matrix with values of mutual information for a randomized
#'        expression matrix where each column is a linearized version of a random
#'        expression matrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
nullLagDist <- function(matrixData,numRep,randomizations, stepSize){
  numGenes <- ncol(matrixData)
  distMat <- matrix(ncol=(numGenes^2),nrow=randomizations)
  for (counts in 1:randomizations) {
    randomMatrixData <- matrix(ncol=numGenes,nrow=nrow(matrixData))
    for (i in 1:ncol(matrixData)){
      randomMatrixData[,i] <- sample(matrixData[,i])
    }
    randomMIMat <- lagMIMat(randomMatrixData,numRep, stepSize)
    randomMIVec <- matrix(randomMIMat,ncol=(numGenes^2))
    distMat[counts,] <- randomMIVec
  }
  return(distMat)
}

optLagDist <- function(matrixData, numRep, stepSize, randomizations, maxStep){
  numGenes <- ncol(matrixData)
  varMat <- matrix(0, nrow = numGenes, ncol = numGenes)
  stepMat <- matrix(nrow = numGenes, ncol = numGenes)
  miMat <- matrix(nrow = numGenes, ncol = numGenes)
  pValMat <- matrix(nrow = numGenes, ncol = numGenes)
  for (i in 1:maxStep){
    step <- stepSize * i
    tempLagMat <- lagMIMat(matrixData, numRep, step)
    tempNullDist <- nullLagDist(matrixData, numRep, randomizations, step)
    tempPMat <- infoScore(tempLagMat, tempNullDist, numGenes, randomizations)
    for (j in 1:numGenes){
      for (k in 1:numGenes){
          if (j>k){
            if (tempPMat[j,k]<=0.05 | tempPMat[k,j]<=0.05){
              locVar <- var(c(tempLagMat[j,k],tempLagMat[j,k]))
              if (locVar > varMat[j,k]){
                varMat[j,k] <- locVar
                stepMat[j,k] <- i
                miMat[j,k] <- tempLagMat[j,k]
                miMat[k,j] <- tempLagMat[k,j]
                pValMat[j,k] <- tempPMat[j,k]
                pValMat[k,j] <- tempPMat[k,j]
              }
            } else {
              miMat[j,k] <- tempLagMat[j,k]
              miMat[k,j] <- tempLagMat[k,j]
              pValMat[j,k] <- tempPMat[j,k]
              pValMat[k,j] <- tempPMat[k,j]
            }
          }
        }
      }
  }
  upperTriangle(varMat)<- lowerTriangle(varMat)
  upperTriangle(stepMat)<- lowerTriangle(stepMat)
  retList <- list(miMat = miMat, pVals = pValMat, stepMat = stepMat)
  return(retList)
}
# 1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables #
option_list = list(
  make_option(c("-e", "--expression_matrix"), type = "character", default = NULL,
              help ="The expression matrix file.", metavar = "character"),
  make_option(c("-i", "--num_reps"), type = "integer", default = 1,
              help ="Number of repetitions per time point. [default= %default]", metavar = "character"),
  make_option(c("-r", "--num_rand"), type = "integer", default = 1000,
              help ="Number of randomizations to perform in order to estimate empirical p values. [default= %default]", metavar = "character"),
  make_option(c("-d", "--dynamical_step"), type = "logical", default = FALSE,
              help ="Perform dynamical step optimization to achieve the optimal step size for each gene pair. [default= %default]", metavar = "character"),
  make_option(c("-m", "--max_step"), type = "numeric", default = 3,
              help ="Maximum step to consider when performing optimization. [default= %default]", metavar = "character"),
  make_option(c("-s", "--step_size"), type = "numeric", default = 1,
              help ="Step size to be lagged (invalid if -d TRUE). [default= %default]", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "edge_list.txt",
              help ="An edge for the esrtimated interactions, includes gene pairs mutual information values and p values  [default= %default]", metavar = "character")
);
# Add command line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt_names <- names(opt)
if (is.null(opt$expression_matrix)){
  print_help(opt_parser)
  err_str <- 'Argument missing "--expression_matrix" must be provided.\n'
  stop(err_str, call.=FALSE)
}
# Parse the command line arguments
expressionMatrix <- opt$expression_matrix
numRep <- opt$num_reps
randomizations <- opt$num_rand
dynamical <- opt$dynamical_step
maxStep <- opt$max_step
stepSize <- opt$step_size
output <- opt$output
# 1.3.2 Load data #
cat('Loading data... ')
matrixData <- read.table(expression_matrix, h = T, stringsAsFactors = F)
numGenes <- ncol(matrixData)
geneNames <- colnames(matrixData)
cat('Done!\n')

#====== 2.0 Estimate pairs of mutual information and their significance ======#
# 2.1 Estimate data sphericity ===============================================#
mlmfit <- lm(as.matrix(iniData)~1)
stopifnot(mauchly.test(mlmfit)$p.value < 0.05)
# 2.2 Calculate mutual information values ====================================#
if (dynamical){
  #If dynamical optimization is specified
  optMutual <- optLagDist(matrixData, numRep, stepSize, randomizations, maxStep)
  adjMI <- optMutual$miMat
  adjPvals <- optMutual$pVals
  adjStep <- optMutual$stepMat
  edgeList <- data.frame(gene1 = c(), gene2 = c(), MI = c(), pVal = c(), step = c())
  for (i in 1:numGenes){
    for (j in 1:numGenes)
      if(i>j){
        locDF <- data.frame (gene1 = geneNames[i], gene2 = geneNames[j], MI = adjMI[i,j], pVal = adjPvals[i,j], step = adjStep[i,j])
        edgeList <- rbind(edgeList, locDF)
        locDF <- data.frame (gene1 = geneNames[j], gene2 = geneNames[i], MI = adjMI[j,i], pVal = adjPvals[j,i], step = adjStep[j,i])
        edgeList <- rbind(edgeList, locDF)
      }
  }
} else {
  iniMutual <- lagMIMat(matrixData,numRep,stepSize)
  mutualNull <- nullLagDist(matrixData,numGenes,randomizations)
  pValues <- infoScore(iniMutual,mutualNull,numGenes,randomizations)
  rownames(iniMutual) <- geneNames
  colnames(iniMutual) <- geneNames
  edgeList <- data.frame(gene1 = c(), gene2 = c(), MI = c(), pVal = c())
  for (i in 1:numGenes){
    for (j in 1:numGenes)
      if(i>j){
        locDF <- data.frame (gene1 = geneNames[i], gene2 = geneNames[j], MI = iniMutual[i,j], pVal = pValues[i,j])
        edgeList <- rbind(edgeList, locDF)
        locDF <- data.frame (gene1 = geneNames[j], gene2 = geneNames[i], MI = iniMutual[j,i], pVal = pValues[j,i])
        edgeList <- rbind(edgeList, locDF)
      }
  }
}

#======================= 3.0 Write the edge list file ========================#
write.table(x = edgeList, file = output, quote = F, sep = '\t', row.names = F)
#=============================================================================#
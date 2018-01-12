mutualInfoEst <- function(matrixData,numGenes) {
    
    numBins <- ceiling(log2(nrow(matrixData)+1))
  mutualMat <- matrix(ncol=numGenes,nrow=numGenes)
    
    for (i in 1:numGenes) {
          for (j in 1:numGenes) {
                  discretVec <- discretize2d(matrixData[,i],matrixData[,j],numBins,numBins)
        mutualMat[i,j] <- mi.empirical(discretVec,unit=c("log2"))
            }
    }
    
    return(mutualMat)
}

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

lagMIMat <- function(matrixData,numRep){
    numGenes <- ncol(matrixData)
  numBins <- ceiling(log2(nrow(matrixData)+1))
    headMatrix <-matrixData[1:(nrow(matrixData)-numRep),]
    tailMatrix <- matrixData[(numRep+1):(nrow(matrixData)),]
      lagMatrix <- matrix(ncol=numGenes,nrow=numGenes)
      for (i in 1:numGenes){
            for (j in 1:numGenes){
                    lagVec<- discretize2d(headMatrix[,i],tailMatrix[,j],numBins,numBins)
            lagMatrix[i,j] <- mi.empirical(lagVec,unit=c("log2"))
                }
        }
        return(lagMatrix)
}

nullLagDist <- function(matrixData,numRep,randomizations){
    numGenes <- ncol(matrixData)
  distMat <- matrix(ncol=(numGenes^2),nrow=randomizations)
    
    for (counts in 1:randomizations) {
          randomMatrixData <- matrix(ncol=numGenes,nrow=nrow(matrixData))
      for (i in 1:ncol(matrixData)){
              randomMatrixData[,i] <- sample(matrixData[,i])
          }
          randomMIMat <- lagMIMat(randomMatrixData,numRep)
          randomMIVec <- matrix(randomMIMat,ncol=(numGenes^2))
              distMat[counts,] <- randomMIVec
            }
    return(distMat)
}


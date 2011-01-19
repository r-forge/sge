###########################################################################/**
# @RdocFunction matrixFromRegions
#
# @title "Data matrix from regions"
#
# 
# \description{
#  @get "title". method to create a matrix for all the genes/snps from 
#   segmented regions
# }
# 
# @synopsis
#
# \arguments{
#   \item{cbsData}{CBSmodel from aroma.affymetrix of GE or CN data.}
#   \item{gi}{Genome information file related to the cbsData.}
# }
#
# \value{
#  \item{dataMat}{LOCIxSAMPLES matrix created from the initial given regions.}
#  }
# \details{
#   The algorithm creates a matrix with the number of rows given by the 
#   genome information file, and the data from the regions given by the 
#   CBSmodel.
# }
# \seealso{
#   Check aroma.affymetrix copy number analysis to better understand what is 
#   being done here.
#   @see "sGEmethod".
# }
#
#*/########################################################################### 

matrixFromRegions <- function(cbsData, gi){
  reg <- cbsData;
  nSamples <- length(reg);
  nbrLoci <- length(getUnitsOnChromosome(gi, c(1:23)));
  dataMat <- matrix(data=0, nrow = nbrLoci, ncol = nSamples);

  nChr <- c(1:23);
  if (nChr[1]==""){
    nChr <- nChr[2:length(nChr)];
  }
  for (i in 1:nSamples){
    sData <- reg[[i]];
    lociDone <- 0;
    for (j in 1:length(nChr)){
      lociChr <- getUnitsOnChromosome(gi, j);
      if  (length(lociChr)==0){
        break;         
      }
      lociPos <- getPositions(gi,units = lociChr);
      lociPos <- sort(lociPos);
      indChrLoci <- sData$chromosome==j;
      scData <- sData[indChrLoci,];
      indReg <- 1;
      nLoci <- length(lociPos);
      auxMat <- matrix(data=0,nrow=nLoci,ncol=1);
      for (k in 1:nLoci){
        locusPos <- lociPos[k];
        found = FALSE;
        while(indReg < sum(indChrLoci) & found == FALSE){
          if(locusPos < scData$start[indReg] & (indReg ==1 || (indReg>1 & locusPos > scData$stop[indReg-1]))){
            auxMat[k,1] <- scData$mean[indReg];
            found = TRUE;    
          }else{
            if(locusPos >= scData$start[indReg] & locusPos <= scData$stop[indReg]){
              auxMat[k,1] <- scData$mean[indReg];
              found = TRUE;
            }else{
              indReg <- indReg + 1;
            }
          }
        }
        if (found == FALSE & indReg == sum(indChrLoci)){
          auxMat[k,1] <- scData$mean[indReg];
        }
      }
      dataMat[(lociDone+1):(lociDone+nLoci),i]  <- auxMat;                   
      lociDone <- lociDone+nLoci;
    }
  }
  return(dataMat);
}
############################################################################
# HISTORY:
# 2011-01-18 [MO]
# o Created.
############################################################################

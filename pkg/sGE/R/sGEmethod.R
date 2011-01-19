###########################################################################/**
# @RdocFunction sGEmethod
#
# @title "Segmentation of GE data"
#
# \description{
#  @get "title". This function does the segmentation of both CN and GE data.
# }
# \arguments{
#   \item{projectName}{String with the name of the project to analyze.}
#   \item{chipTypes}{List with the names of both (GE and CN) cdfs.}
#   \item{fncCN}{Name of the function to use to summarize the CN data.}
# }
# \value{
#  Returns a @list:
#  \item{matCNsnps}{The SNPSxSAMPLES @matrix containing segmented CN estimates.}
#  \item{matCNgenes}{A GENESxSAMPLES @matrix containing segmented CN estimates.}
#  \item{matGEgenes}{A GENESxSAMPLES @matrix containing segmented GE estimates.}
#  }
# \details{
#   The algorithm applies CBS segmentation to both CN and GE data. It also
#   generates a GENESxSAMPLES matrix of CN estimates in order to make it more
#   comparable to GE data.
# }
# \seealso{
#   Check aroma.affymetrix copy number analysis to better understand what is 
#   being done here.
#   @see "matrixFromRegions".
# }
#*/########################################################################### 
sGEmethod <- function(projectName=NULL, chipTypes=NULL, fncCN="avg"){
  if(is.null(projectName) | is.null(chipTypes)){
    throw("not enough arguments for the function");    
  }
  if(length(chipTypes) != 2){
    throw("wrong number of chipTypes");
  }
  if(fncCN != "avg" & fncCN != "rma"){
    throw("wrong summarization method for CN data");
  }
  
  verbose <- Arguments$getVerbose(-8);
  timestampOn(verbose);
  
  ### CN analysis       
  cdfCN <- AffymetrixCdfFile$byChipType(chipTypes[[1]]);
  print(cdfCN);
  
  giCN <- getGenomeInformation(cdfCN);
  print(giCN);
  
  siCN <- getSnpInformation(cdfCN);
  print(siCN);

  csCN <- AffymetrixCelSet$byName(projectName, cdf=cdfCN);
  print(csCN);
  
  accCN <- AllelicCrosstalkCalibration(csCN, model="CRMAv2");
  print(accCN);
  csCcn <- process(accCN, verbose=verbose);

  bpnCN <- BasePositionNormalization(csCcn, target="zero");
  print(bpnCN);
  csNcn <- process(bpnCN, verbose=verbose);

  #using AVG
  if (fcnCN=="avg"){
    plmCN <- AvgCnPlm(csNcn, mergeStrands=TRUE, combineAlleles=TRUE);
  }else{
    plmCN <- RmaCnPlm(csNcn, mergeStrands=TRUE, combineAlleles=TRUE);
  }
  fit(plmCN , verbose=verbose)
  
  # Getting gene expression values
  cesCN <- getChipEffectSet(plmCN)
  # Segmentation
  cbsCN <- CbsModel(cesCN)
  fit(cbsCN ,verbose = verbose)
  regionsCN <- getRegions(cbsCN)
  # Creation of a matrix from the regions
  matCNsnps <- matrixFromRegions(regionsCN, giCN)

  #### GE analysis

  cdfGE <- AffymetrixCdfFile$byChipType(chipTypes[[2]]);
  
  csGE <- AffymetrixCelSet$byName(projectName , cdf=cdfGE )
  giGE <- getGenomeInformation (cdfGE ,verbose =verbose)
  # Background removal
  BCge <- RmaBackgroundCorrection(csGE)
  csBCge <- process(BCge, verbose = verbose )
  # Normalization
  qnGE <- QuantileNormalization(csBCge , typesToUpdate="pm")
  csNge <- process(qnGE , verbose=verbose )
  # Summarization
  plmGE <- RmaPlm(csNge)
  fit(plmGE , verbose=verbose)
  # Getting gene expression value s
  cesGE <- getChipEffectSet(plmGE)
  # Segmentation
  cbsGE <- CbsModel(cesGE)
  fit(cbsGE ,verbose = verbose)
  regionsGE <- getRegions(cbsGE)
  # Creation of a matrix from the regions
  matGEgenes <- matrixFromRegions(regionsGE, giGE)
  matCNgenes <- matrixFromRegions(regionsCN, giGE)  

  return(list(matCNsnps, matCNgenes, matGEgenes));
}
############################################################################
# HISTORY:
# 2011-01-18 [MO]
# o Created.
############################################################################


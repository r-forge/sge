%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  matrixFromRegions.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{matrixFromRegions}
\alias{matrixFromRegions}


 \title{Data matrix from regions}


 \description{
  Data matrix from regions. method to create a matrix for all the genes/snps from
   segmented regions
 }

 \usage{matrixFromRegions(cbsData, gi)}

 \arguments{
   \item{cbsData}{CBSmodel from aroma.affymetrix of GE or CN data.}
   \item{gi}{Genome information file related to the cbsData.}
 }

 \value{
  \item{dataMat}{LOCIxSAMPLES matrix created from the initial given regions.}
  }
 \details{
   The algorithm creates a matrix with the number of rows given by the
   genome information file, and the data from the regions given by the
   CBSmodel.
 }
 \seealso{
   Check aroma.affymetrix copy number analysis to better understand what is
   being done here.
   \code{\link{sGEmethod}}().
 }


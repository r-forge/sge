%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  sGEmethod.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{sGEmethod}
\alias{sGEmethod}


\title{Segmentation of GE data}

\description{
 Segmentation of GE data. This function does the segmentation of both CN and GE data.
}
\arguments{
  \item{projectName}{String with the name of the project to analyze.}
  \item{chipTypes}{List with the names of both (GE and CN) cdfs.}
  \item{fncCN}{Name of the function to use to summarize the CN data.}
}
\value{
 Returns a \code{\link[base]{list}}:
 \item{matCNsnps}{The SNPSxSAMPLES \code{\link[base]{matrix}} containing segmented CN estimates.}
 \item{matCNgenes}{A GENESxSAMPLES \code{\link[base]{matrix}} containing segmented CN estimates.}
 \item{matGEgenes}{A GENESxSAMPLES \code{\link[base]{matrix}} containing segmented GE estimates.}
 }
\details{
  The algorithm applies CBS segmentation to both CN and GE data. It also
  generates a GENESxSAMPLES matrix of CN estimates in order to make it more
  comparable to GE data.
}
\seealso{
  Check aroma.affymetrix copy number analysis to better understand what is
  being done here.
  \code{\link{matrixFromRegions}}().
}


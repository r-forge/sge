library("aroma.affymetrix");
library("sGE");

projectName <- "sGE,test"
chipTypes <- list("Mapping250K_Nsp","HGU133Plus2_Hs_ENSG");

dataList <- sGEmethod(projectName, chipTypes, fncCN="rma");
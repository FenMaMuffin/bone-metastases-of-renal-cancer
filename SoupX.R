library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(SoupX)
library(Rnmr1D)

run_SoupX<-function(dataDir,geneSet){
  sc<-load10X(dataDir)
  plt<-plotMarkerMap(sc,geneSet)+ ggtitle(paste0(geneSet,collapse='_'))
  sc<-autoEstCont(sc)
  sc_adj<-adjustCounts(sc,roundToInt = T)
  seurat_obj<-CreateSeuratObject(sc_adj)
  return(list(plt=plt,seurat_obj=seurat_obj))
}


setwd("~/kidneycancer/SoupX")


Dir='~/bone_metastatic/rawdata/P01_P/outs/'
soupX_out<-run_SoupX(dataDir= Dir,geneSet = 'IGKC')
saveRDS(soupX_out$seurat_obj,"P01_P.rds")
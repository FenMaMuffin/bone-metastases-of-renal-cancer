#doubletFinder
library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(SoupX)
library(Rnmr1D)

run_doubletfinder<-function(object,n_capture_cell){
  n_capture_cell<-as.character(n_capture_cell)
  multiplet_rate<-c("500"=0.004,'1000'=0.008,'2000'=0.016,'3000'=0.023,'4000'=0.031,
                    '5000'=0.039,'6000'=0.046,'7000'=0.054,'8000'=0.061,'9000'=0.069,
                    '10000'=0.076)
  stopifnot(n_capture_cell%in%names(multiplet_rate))
  
  ##pre-process Seurat object
  seurat_obj<-SCTransform(object)%>%RunPCA%>%RunUMAP(dim=1:30)%>%FindNeighbors(reduction = "pca", dims = 1:30)%>%FindClusters(resolution = 0.8)
  
  ##pK identification
  n_cpu<-round(parallel::detectCores()*0.5)
  sweep.res.list<-paramSweep_v3(seurat_obj,PCs = 1:30,sct = T,num.cores = n_cpu)
  sweep.stats<-summarizeSweep(sweep.res.list,GT=F)
  bcmvn<-find.pK(sweep.stats)
  
  param<-c(pN=0.25,pK=bcmvn$pK[which.max(bcmvn$BCmetric)]%>%as.character()%>%as.numeric(),
           nExp=round(multiplet_rate[n_capture_cell]*nrow(seurat_obj@meta.data))%>%unname())
  #nExp=round(round(multiplet_rate[n_capture_cell]*ncol(seurat_obj)) * (1-modelHomotypic(seurat_obj$seurat_clusters)))
  
  seurat_obj<-doubletFinder_v3(seurat_obj,PCs=1:30,pN=param['pN'],pK=param['pK'],nExp = param['nExp'],
                               reuse.pANN = F,sct=T)
  seurat_obj$Is_doublet<-seurat_obj[[paste0('DF.classifications_',param['pN'],'_',param['pK'],'_',param['nExp'])]]=='Doublet'
  
  plt<-DimPlot(seurat_obj,group.by = paste0('DF.classifications_',param['pN'],'_',param['pK'],'_',param['nExp']))
  
  
  return(list(plt=plt,seurat_obj=seurat_obj))
}


setwd("~/kidneycancer/DoubletFinder")

datalist<-c("P01_P.rds",  "P02_P.rds",  "P03_P.rds", "P04_P.rds", "P05_P.rds", "P06_P.rds")
for (i in 1:length(datalist)){
  current_dir='~/kidneycancer/01.data/primaryccRCC/'
  input_dir=paste0(current_dir,datalist[i])
  soupx_out<-readRDS(input_dir)
  if(ncol(soupx_out)<1000){
    doubletfinder_out<-run_doubletfinder(object = soupx_out,n_capture_cell = 500)}
  else{doubletfinder_out<-run_doubletfinder(object = soupx_out,n_capture_cell = 5000)}
  saveRDS(doubletfinder_out,datalist[i])
}
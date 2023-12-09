library(Seurat)
library(dplyr)
library(harmony)

setwd('~/kidneycancer/Output')
current_dir='~/kidneycancer/processed_data/'
data<-readRDS(paste0(current_dir,"Mergedata_raw.rds"))
data<-subset(data,DoubletFinder=='Singlet')

mtGenes=grep('^MT-',rownames(data),value=TRUE)
red = data[!rownames(data) %in% mtGenes,]
mtFrac = colSums(data[mtGenes,])/colSums(data)
red@meta.data$percent.MT = mtFrac


HSPGenes=grep('^HSP',rownames(data),value=TRUE)
red2 = red[!rownames(red) %in% HSPGenes,]
HSPFrac = colSums(data[HSPGenes,])/colSums(data)
red2@meta.data$percent.HSP = HSPFrac

RPGenes=grep('^RP[SL]',rownames(data),value=TRUE)
red3 = red2[!rownames(red2) %in% RPGenes,]
RPFrac = colSums(data[RPGenes,])/colSums(data)
red3@meta.data$percent.RP = RPFrac

pdf("Vlnplot_before_1.pdf",width = 15,height = 6)
plot1<-VlnPlot(object = red3,group.by = 'orig.ident',features= c("percent.MT"),pt.size = 0)+NoLegend()
plot2<-VlnPlot(object = red3,group.by = 'orig.ident',features= c("nFeature_RNA"),pt.size = 0)+NoLegend()
plot3<-VlnPlot(object = red3,group.by = 'orig.ident',features= c("nCount_RNA"),pt.size = 0)+NoLegend()
CombinePlots(plots = list(plot1,plot2,plot3))
dev.off()
pdf("Vlnplot_before_2.pdf",width = 15,height = 3)
plot1<-VlnPlot(object = red3,group.by = 'orig.ident',features= c("percent.HSP"),pt.size = 0)+NoLegend()
plot2<-VlnPlot(object = red3,group.by = 'orig.ident',features= c("percent.RP"),pt.size = 0)+NoLegend()
CombinePlots(plots = list(plot1,plot2))
dev.off()


#SCTransform
seurat<-subset(red3,subset=nFeature_RNA> 700 & nFeature_RNA< 8000 & nCount_RNA<20000 & percent.MT<0.1)
seurat<- SCTransform(seurat, vars.to.regress=c("percent.MT",'percent.HSP','percent.RP'),verbose = FALSE)

seurat<- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:25, verbose = FALSE)
seurat<- FindNeighbors(seurat, dims = 1:25, verbose = FALSE)
seurat <- FindClusters(seurat, verbose = FALSE)

dir="~/kidneycancer/processed_data/"
saveRDS(seurat,paste0(dir,'data_SCT.rds'))

dir='~/kidneycancer/Output/'
pdf(paste0(dir,"Vlnplot_after_1.pdf"),width = 15,height = 6)
plot1<-VlnPlot(object = seurat,group.by = 'orig.ident',features= c("percent.MT"),pt.size = 0)+NoLegend()
plot2<-VlnPlot(object = seurat,group.by = 'orig.ident',features= c("nFeature_RNA"),pt.size = 0)+NoLegend()
plot3<-VlnPlot(object = seurat,group.by = 'orig.ident',features= c("nCount_RNA"),pt.size = 0)+NoLegend()
CombinePlots(plots = list(plot1,plot2,plot3))
dev.off()
pdf(paste(dir,"Vlnplot_after_2.pdf"),width = 15,height = 3)
plot1<-VlnPlot(object = seurat,group.by = 'orig.ident',features= c("percent.HSP"),pt.size = 0)+NoLegend()
plot2<-VlnPlot(object = seurat,group.by = 'orig.ident',features= c("percent.RP"),pt.size = 0)+NoLegend()
CombinePlots(plots = list(plot1,plot2))
dev.off()
plot<-ElbowPlot(seurat, ndims = 50)
pdf(paste0(dir,'Elbowplot_SCT.pdf'),width = 6,height = 6)
plot
dev.off()
pdf(paste0(dir,"UMAP_SCT.pdf"),width =6,height= 6)
DimPlot(seurat, label = TRUE) + NoLegend()
dev.off()
pdf(paste0(dir,"groupcluster_SCT_1.pdf"),width =6,height= 6)
DimPlot(seurat, label = TRUE, group.by= "orig.ident") + NoLegend()
dev.off()
pdf(paste0(dir,"groupcluster_SCT_2.pdf"),width =6,height= 6)
DimPlot(seurat, label = TRUE, group.by= "Group") + NoLegend()
dev.off()
pdf(paste0(dir,"PCA_SCT_1.pdf"),width=8,height = 6)
DimPlot(seurat,dims = 1:2,reduction = "pca",group.by= "orig.ident")
dev.off()
pdf(paste0(dir,"PCA_SCT_2.pdf"),width=8,height = 6)
DimPlot(seurat,dims = 1:2,reduction = "pca",group.by= "Group")
dev.off()


#Harmony
seurat<-RunHarmony(seurat,group.by.vars='orig.ident',assay.use="SCT",max.iter.harmony=20,lambda=1)
seurat<-RunUMAP(seurat,reduction="harmony", dims=1:30)
seurat<-FindNeighbors(seurat, dims = 1:25) %>% FindClusters()

dir="~/kidneycancer/processed_data/"
saveRDS(seurat,paste0(dir,'data_SCT_harmony.rds'))

pdf("UMAP_SCT_Harmony.pdf",width =6,height= 6)
DimPlot(seurat, label = TRUE) + NoLegend()
dev.off()

logFCfilter=0.5  
adjPvalFilter=0.05  
markers<-FindAllMarkers(object=seurat,only.pos=FALSE,min.pct=0.25,logfc.threshold=logFCfilter) 
sig.marker=markers[abs(as.numeric(as.vector(markers$avg_logFC)))>logFCfilter&as.numeric(as.vector(markers$p_val_adj))<adjPvalFilter]
write.csv(pbmc.markers,'Markergene_SCT_harmony.csv')

top10<-markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_logFC)  
pdf(file="markerHeatmap_SCT_harmony.pdf",width=20,height=9)
DoHeatmap(object=seurat,features=top10$gene)+NoLegend()
dev.off()

markergene<-c("CD3D","CD3E","GNLY","FGFBP2","CD14","FCGR3A","LYZ","XCL2","KLRC1","PLVAP","PECAM1","KRT18","VEGFA","COL1A1","COL1A2","CD79A","MS4A1","TPSB2","TPSAB1","CXCL8","S100A8","S100A9")
pdf("Markergene_SCT_Harmony.pdf",width=16,height=12)
FeaturePlot(seurat,features=markergene,cols = c("#B4CDCD", "#FF6A6A"))
dev.off()
pdf('markergene_dotplot_SCT_harmony.pdf')
DotPlot(seurat,features = markergene)
dev.off()

name<-read.csv(paste0(dir,'Rename_SCT_harmony.csv'),header = T)
new_iden<-name$Anntation
names(new_iden)<- levels(seurat)
seurat2<-RenameIdents(seurat,new_iden) 
seurat2@meta.data$name<-seurat2@active.ident
saveRDS(seurat2,paste0(dir,'data_SCT_harmony_rename.rds'))
out_dir='~/kidneycancer/Output/'
pdf(paste0(out_dir,'Umap_SCT_harmony_rename.pdf'),width = 8,height = 8)
DimPlot(seurat2, label = TRUE) + NoLegend()
dev.off()


dir='~/kidneycancer/processed_data/'
alldata<-readRDS(paste0(dir,'data_SCT_harmony_rename.rds'))

celltype<-c("TNK cells",'Myeloid cells','Cancer cells','Fibroblasts','Endothelial cells','B cells','Mast cells',"Plasma cells")

out_dir='~/kidneycancer/processed_data/SCT_harmony/'
for (i in 1:length(celltype)) {
  seurat<-subset(alldata,idents=celltype[i])
  DefaultAssay(seurat)<-'RNA'
  seurat <- NormalizeData(seurat) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData(vars.to.regress=c('percent.MT','percent.HSP','percent.RP'))
  all.genes <- rownames(seurat)
  seurat <- ScaleData(seurat, features = all.genes)
  seurat <- RunPCA(seurat, verbose = FALSE)  
  seurat<-RunHarmony(seurat,group.by.vars='orig.ident',assay.use="RNA",max.iter.harmony=20,lambda=1)
  seurat<-RunUMAP(seurat,reduction="harmony", dims=1:25)
  seurat<-FindNeighbors(seurat, dims = 1:25) %>% FindClusters()
  saveRDS(seurat,paste0(out_dir,celltype[i],'.rds'))
}




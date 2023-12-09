library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)

setwd("~/Renal cancer/Myleid/monocle")
myeloid<-readRDS("myeloid_final.rds")
meta <- myeloid@meta.data
meta$active.ident<-myeloid@active.ident
exp <- myeloid@assays$RNA@counts
genes <- as.data.frame(rownames(exp))
colnames(genes)<-"gene_short_name"
rownames(genes) <- genes$gene_short_name

set.seed(6)
downsample_cells <- as.vector(sample_n(as.data.frame(myeloid), 2000,replace = FALSE)$Macro)
exp_data_down <- exp[,downsample_cells]
meta_data_down <- meta[downsample_cells,]

pd <- new("AnnotatedDataFrame", data = meta_data_down)	
fd <- new("AnnotatedDataFrame", data = genes)
HSMM <- newCellDataSet(exp_data_down,
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1) 

clustering_DEG_genes <- differentialGeneTest(HSMM, fullModelFormulaStr = '~active.ident', cores = 8)
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:400]
HSMM1<-HSMM
HSMM <- setOrderingFilter(HSMM, ordering_genes = HSMM_ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

plot_cell_trajectory(HSMM, color_by = "active.ident",cell_size=1.5)+theme(legend.position = "right")+scale_color_manual(values = colors)+guides(color = guide_legend(override.aes = list(size=5)))
pData(HSMM)$Pseudotime <- max(pData(HSMM)$Pseudotime) - pData(HSMM)$Pseudotime
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
df <- pData(HSMM)
ggplot(df, aes(Pseudotime, colour = active.ident, fill=active.ident)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+ facet_wrap(~active.ident, nrow = 4)+scale_color_manual(values = colors)
plot_pseudotime_heatmap(HSMM[HSMM_ordering_genes,], num_clusters = 3, cores = 1, show_rownames = T)



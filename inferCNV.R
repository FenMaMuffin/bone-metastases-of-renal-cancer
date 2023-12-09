library(Seurat)
library(ggplot2)
library(infercnv)

setwd('~/kidneycancer/IntergateOutput/Cancer/inferCNV/count')

expFile="expFile_count_group.txt"
groupFile="groupFile_group.txt"
geneFile="geneFile.txt"
infercnv_obj<-CreateInfercnvObject(raw_counts_matrix=expFile,
                                   annotations_file=groupFile,
                                   delim="\t",
                                   gene_order_file=geneFile,                               ref_group_names=c("Fibroblasts",'Endothelial cells'))
infercnv_all<-infercnv::run(infercnv_obj,
                            cutoff=0.1,
                            out_dir="~/kidneycancer/IntergateOutput/Cancer/inferCNV/count/CNVoutput",
                            cluster_by_groups=T,
                            num_threads=30,
                            denoise=F,
                            HMM=F )

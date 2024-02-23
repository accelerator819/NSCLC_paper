library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(clusterProfiler)
setwd('/mnt/d/project_ptc/PTC-lung/single_cell/scanpy/paper_scripts')
#Convert("./scanpy_rawcount_7734_0129.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat <- LoadH5Seurat("./scanpy_rawcount_7734_0129.h5seurat")
annodata <- read.csv('./merged_cellAnnotation.csv',row.names = 1)
annodata <- annodata[row.names(seurat@meta.data),]
seurat@meta.data <- annodata
#row.names(adata@meta.data) == row.names(annodata)
#View(adata@meta.data)
#View(annodata)
#rm(adata)
#rawcount <- read.csv('./rawcounts_merged7734_0129.csv',row.names = 1)

counts <- seurat@assays$RNA@counts 

metadata <- seurat@meta.data
seurat@active.ident <- as.factor(seurat$major_clusters)
# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)
seurat$sampleId <- str_replace_all(seurat$sampleId, c("_" = ""))
metadata$sample_id <- factor(seurat$sampleId)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

# Identify groups for aggregation of counts
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Named vector of cluster names
kids <- purrr::set_names(levels(sce$cluster_id))
kids

# Total number of clusters
nk <- length(kids)
nk

# Named vector of sample names
sids <- purrr::set_names(levels(sce$sample_id))
sids
# Total number of samples 
ns <- length(sids)
ns

# Generate sample level metadata

## Determine the number of cells per sample
table(sce$sample_id)

## Turn named vector into a numeric vector of number of cells per sample
n_cells <- as.numeric(table(sce$sample_id))

## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
m <- match(sids, sce$sample_id)

## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"cluster_id")
ei <- ei[,c('sampleType','sample_id','n_cells')]

# Aggregate the counts per sample_id and cluster_id

# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "sample_id")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 
pb <- t(pb)
pbulk <-  as.matrix(pb)[,c(-17,-18)]
cn <- str_replace_all(colnames(pbulk),'p0129ori','0129-T')
cn <- str_replace_all(cn,'p0129ptc','0129-C')
cn <- str_replace_all(cn,'p7734ori','7734-T')
cn <- str_replace_all(cn,'p7734ptc','7734-C')
colnames(pbulk) <- cn
coldata <- data.frame(row.names = colnames(pbulk),
                      celltype = c(rep('dc',4),rep('epi',4),
                                   rep('fib',4),rep('mo',4)))
dds <- DESeqDataSetFromMatrix(pbulk, 
                              colData = coldata, 
                              design = ~ celltype)
dds <- DESeq(dds)
ncount <- counts(dds,normalized = T)
ncount_0129 <- ncount[,c(1,2,5,6,9,10,13,14)]

ncount_7734 <- ncount[,c(3,4,7,8,11,12,15,16)]
cor_mtx_0129 <-cor(log2(ncount_0129+1),method = "spearman")
dir <- '/mnt/d/project_ptc/PTC-lung/single_cell/scanpy/paper_scripts'
png(file.path(dir,'0129_scRNACorrelation.png'),res = 600,height = 8,width = 8,units = 'in')
pheatmap(cor_mtx_0129,cellwidth = 30,cellheight = 30,
         cluster_rows = F,cluster_cols = F) 
dev.off()

cor_mtx_7734 <-cor(log2(ncount_7734+1),method = "spearman")
dir <- '/mnt/d/project_ptc/PTC-lung/single_cell/scanpy/paper_scripts'
png(file.path(dir,'7734_scRNACorrelation.png'),res = 600,height = 8,width = 8,units = 'in')
pheatmap(cor_mtx_7734,cellwidth = 30,cellheight = 30,
         cluster_rows = F,cluster_cols = F) 
dev.off()



library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(tidyverse)
library(tximport)
library(sva)
library(RColorBrewer)
### Scripts description ###
### This scripts is intended for integrating ptc-tissue RNA-expression data with TCGA###
### and peform correlation analysis ###
dir = 'd:/project_ptc/PTC-lung/RNAexp'
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
workflow_type <- "HTSeq - Counts"
query_TranscriptomeCounts <- GDCquery(project = "TCGA-LUAD",
                                      data.category = data_category,
                                      data.type = data_type,
                                      workflow.type = workflow_type)
#GDCdownload(query_TranscriptomeCounts,method = 'api')
luad_exp = GDCprepare(query = query_TranscriptomeCounts,
                      directory = 'd:/project_ptc/PTC-lung/GDCdata')
luad_exp = luad_exp[,is.na(luad_exp$paper_expression_subtype) == FALSE]

query_TranscriptomeCounts <- GDCquery(project = "TCGA-LUSC",
                                      data.category = data_category,
                                      data.type = data_type,
                                      workflow.type = workflow_type)
#GDCdownload(query_TranscriptomeCounts,method = 'api')
lusc_exp = GDCprepare(query = query_TranscriptomeCounts,
                      directory = 'd:/project_ptc/PTC-lung/GDCdata')
lusc_exp = lusc_exp[,is.na(lusc_exp$paper_Expression.Subtype) == FALSE]





#### PTC-tissue rawcount
sample_info = read.table(file.path(dir,'sample_meta.txt'),sep = '\t',header = T)
row.names(sample_info) = sample_info$sample_name
sample_info = sample_info[sample_info$condition!='tissue',]
files = file.path(dir,sample_info$file_name,"quant.sf")
names(files) = sample_info$sample_name
tx2gene = read.csv('d:/project_ptc/tx2gene_v96.csv',header = T)
txi = tximport(files,type="salmon",tx2gene=tx2gene[,c(1,2)],ignoreTxVersion = TRUE,countsFromAbundance="lengthScaledTPM")
raw_counts = round(txi$counts)
colnames(raw_counts) = sample_info$sample_name
raw_counts = as.data.frame(raw_counts)%>% add_column(gene_id=row.names(raw_counts),.before = 1)
anno_table = dplyr::filter(tx2gene, gene_id %in% row.names(raw_counts))
anno_table = anno_table[which(duplicated(anno_table$symbol) == FALSE),]
raw_counts = merge(anno_table,raw_counts,by.x = 'gene_id',by.y = 'gene_id')
row.names(raw_counts) = raw_counts$symbol
raw_counts = raw_counts[,c(-1,-2,-3,-4,-5)]
meta_custom = data.frame(row.names = sample_info$sample_name,
                         subtype = 'notAvailable',
                         smoking.status = 'notAvailable',
                         cancertype = sample_info$class,
                         batch = 'PTC-TISSUE'
)
#### TCGA-LUAD
luad_count = assay(luad_exp)
luad_count = as.data.frame(luad_count)%>% add_column(gene_id=row.names(luad_count),.before = 1)
anno_table = dplyr::filter(tx2gene, gene_id %in% row.names(luad_count))
anno_table = anno_table[which(duplicated(anno_table$symbol) == FALSE),]
luad_count = merge(anno_table,luad_count,by.x = 'gene_id',by.y = 'gene_id')
row.names(luad_count) = luad_count$symbol
luad_count = luad_count[,c(-1,-2,-3,-4,-5)]
meta_luad = data.frame(row.names = luad_exp$barcode,
                       subtype = luad_exp$paper_expression_subtype,
                       smoking.status = luad_exp$paper_Smoking.Status,
                       cancertype = 'LUAD.TCGA',
                       batch = 'TCGA'
)

#### TCGA-LUSC
lusc_count = assay(lusc_exp)
lusc_count = as.data.frame(lusc_count)%>% add_column(gene_id=row.names(lusc_count),.before = 1)
anno_table = dplyr::filter(tx2gene, gene_id %in% row.names(lusc_count))
anno_table = anno_table[which(duplicated(anno_table$symbol) == FALSE),]
lusc_count = merge(anno_table,lusc_count,by.x = 'gene_id',by.y = 'gene_id')
row.names(lusc_count) = lusc_count$symbol
lusc_count = lusc_count[,c(-1,-2,-3,-4,-5)]
meta_lusc = data.frame(row.names = lusc_exp$barcode,
                       subtype = lusc_exp$paper_Expression.Subtype,
                       smoking.status = lusc_exp$paper_Smoking.Status,
                       cancertype = 'LUSC.TCGA',
                       batch = 'TCGA'
                       )


meta_merged = rbind(meta_custom,meta_luad,meta_lusc)

 #### library size normalization ####

lib.size <- estimateSizeFactorsForMatrix(raw_counts)
normalized_custom <- t(t(raw_counts)/lib.size)
normalized_custom <- as.data.frame(normalized_custom)%>% add_column(gene_id=row.names(normalized_custom),.before = 1)
lib.size <- estimateSizeFactorsForMatrix(luad_count)
normalized_luad <- t(t(luad_count)/lib.size)
normalized_luad <- as.data.frame(normalized_luad)%>% add_column(gene_id=row.names(normalized_luad),.before = 1)
lib.size <- estimateSizeFactorsForMatrix(lusc_count)
normalized_lusc <- t(t(lusc_count)/lib.size)
normalized_lusc <- as.data.frame(normalized_lusc)%>% add_column(gene_id=row.names(normalized_lusc),.before = 1)

df1 = merge(normalized_custom,normalized_luad,by.x = 'gene_id',by.y = 'gene_id')
normalized_merged = merge(df1,normalized_lusc,by.x = 'gene_id',by.y = 'gene_id')
rm(normalized_custom,normalized_luad,normalized_lusc,df1)
rm(raw_counts,lusc_count,luad_count)
row.names(normalized_merged) = normalized_merged$gene_id
normalized_merged = normalized_merged[,-1]
#### correlation

batch = meta_merged$batch
adjusted <- ComBat(as.matrix(log2(normalized_merged+1)),batch = batch,par.prior = T,prior.plots = TRUE)
data.var <- apply(adjusted, 1, stats::var)
most.var.normalized <- adjusted[order(data.var, decreasing = TRUE)[1:500],]
write.table(row.names(adjusted[order(data.var, decreasing = TRUE)[1:500],]),file.path(dir,'500_mostvar.txt'),sep = '\t',row.names = F)


heatmap_annotation = data.frame(row.names = row.names(meta_merged),
                                Data.Source = meta_merged$cancertype)
ann_colors = list(Data.Source = c(LUAD.TCGA = "#ffe4c4",LUSC.TCGA = "#d3eaec",
                                  LUAD.PTC = "#a06d3a",LUSC.PTC = "#2f495e"))
png(file.path(dir,'TCGA-PTC.onlyPTC.png'),height = 10000,width = 10000,res = 600)
pheatmap(most.var.normalized,
         annotation_col = heatmap_annotation,
         annotation_colors = ann_colors,
         show_rownames = F, show_colnames = F,
         colorRampPalette(rev(brewer.pal(n = 10, name ="RdBu")))(100),
         fontsize_col = 2)
dev.off()



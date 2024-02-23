library(tximport)
library(pheatmap)
tx2gene <- read.csv('d:/project_ptc/tx2gene_v96.csv',header = T)
dir <- '/mnt/d/project_ptc/PTC-lung/cluster'
sample_info <- read.table(file.path(dir,"sample_info.txt"), header=TRUE ,sep = '\t')
files <- file.path(dir,sample_info$file_name,'salmon_quant',sample_info$sample_id,"quant.sf")
names(files) <- sample_info$tpm_id
txi = tximport(files,type="salmon",tx2gene=tx2gene[,c(1,2)],ignoreTxVersion = TRUE,countsFromAbundance="lengthScaledTPM")

tpm_counts <- txi$abundance
colnames(tpm_counts) = sample_info$tpm_id
tpm_counts = as.data.frame(tpm_counts) %>% add_column(gene_id = row.names(tpm_counts),.before = 1)
anno_table = dplyr::filter(tx2gene, gene_id %in% row.names(tpm_counts))
anno_table = anno_table[which(duplicated(anno_table$symbol) == FALSE),]
tpm_counts = merge(anno_table,tpm_counts,by.x = 'gene_id',by.y = 'gene_id')
row.names(tpm_counts) = tpm_counts$symbol
tpm_counts = tpm_counts[,c(-1,-2,-3,-4,-5)]
write.table(tpm_counts,file.path(dir,'lung_cluster_tpm.txt'),sep='\t',row.names = T)
genes<-c('EPCAM','CDH1','SFN','CD44','S100A4','PROM1','ALDH1A1','SOX2','FAP','COL6A2',
         'CALD1','PTPRC','CD14','CD68','CD3D')
annotation <- data.frame(row.names = colnames(tpm_counts),
                         source = sample_info$cancer_type_special,
                         patient = sample_info$patient_id)
p <- pheatmap(log2(tpm_counts+1)[genes,],
         annotation_col = annotation,
         show_colnames = F,
         cluster_rows = F, cluster_cols = F)
ggsave(file.path(dir,"marker.png"),p,width=8,height=8,units='in',dpi=400,scale = 1)

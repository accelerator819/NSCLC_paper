dev.off()
rm(list=ls())
library(tidyverse)
library(pheatmap)
dir='/mnt/d/project_ptc/PTC-lung/single_cell/scanpy/paper_scripts/out' 
##########  7734 PE ######
#tumor
mypvals <- read.table(file.path(dir,"7734_tumor","pvalues.txt"), sep='\t',header = T,check.names = F)
mymeans <- read.table(file.path(dir,"7734_tumor","means.txt"), sep='\t',header = T,check.names = F)
goi <- read.table(file.path(dir,'7734_tumor','GOI.txt'),sep='\t')$V1
coi <- read.table(file.path(dir,'7734_tumor','COI.txt'),sep='\t')$V1

means_df_tumor <- mymeans[mymeans$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(means_df_tumor) <- means_df_tumor$interacting_pair
means_df_tumor <- means_df_tumor[,-1]
pval_df_tumor <- mypvals[mypvals$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(pval_df_tumor) <- pval_df_tumor$interacting_pair
pval_df_tumor <- pval_df_tumor[,-1]

label_mtx <- as.matrix(pval_df_tumor)
label_mtx[label_mtx<=0.01] <- '\u2217'
label_mtx[label_mtx>0.01] <- ''
png(file.path(dir,'7734-Pleural Effution.LR.sig.png'),res = 600,height = 7,width = 5,units = 'in')
pheatmap(log2(means_df_tumor+1),cluster_rows = F,cluster_cols = F,display_numbers = label_mtx,
         cellwidth = 25,cellheight = 25,scale='row',main = '7734 Pleural effusion')
dev.off()
# ptc
mypvals <- read.table(file.path(dir,"7734_ptc","pvalues.txt"), sep='\t',header = T,check.names = F)
mymeans <- read.table(file.path(dir,"7734_ptc","means.txt"), sep='\t',header = T,check.names = F)
goi <- read.table(file.path(dir,'7734_tumor','GOI.txt'),sep='\t')$V1
coi <- read.table(file.path(dir,'7734_tumor','COI.txt'),sep='\t')$V1

means_df_ptc <- mymeans[mymeans$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(means_df_ptc) <- means_df_ptc$interacting_pair
means_df_ptc <- means_df_ptc[,-1]
pval_df_ptc <- mypvals[mypvals$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(pval_df_ptc) <- pval_df_ptc$interacting_pair
pval_df_ptc <- pval_df_ptc[,-1]

label_mtx <- as.matrix(pval_df_ptc)
label_mtx[label_mtx<=0.01] <- '\u2217'
label_mtx[label_mtx>0.01] <- ''
png(file.path(dir,'7734-PTC.LR.sig.png'),res = 600,height = 7,width = 5,units = 'in')
pheatmap(log2(means_df_ptc+1),cluster_rows = F,cluster_cols = F,display_numbers = label_mtx,
         cellwidth = 25,cellheight = 25,scale = 'row',main = '7734 PTC')
dev.off()




##########  0129 Tissue #######
#tumor
mypvals <- read.table(file.path(dir,"0129_tumor","pvalues.txt"), sep='\t',header = T,check.names = F)
mymeans <- read.table(file.path(dir,"0129_tumor","means.txt"), sep='\t',header = T,check.names = F)
goi <- read.table(file.path(dir,'0129_tumor','GOI.txt'),sep='\t')$V1
coi <- read.table(file.path(dir,'0129_tumor','COI.txt'),sep='\t')$V1

means_df_tumor <- mymeans[mymeans$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(means_df_tumor) <- means_df_tumor$interacting_pair
means_df_tumor <- means_df_tumor[,-1]
pval_df_tumor <- mypvals[mypvals$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(pval_df_tumor) <- pval_df_tumor$interacting_pair
pval_df_tumor <- pval_df_tumor[,-1]

label_mtx <- as.matrix(pval_df_tumor)
label_mtx[label_mtx<=0.01] <- '\u2217'
label_mtx[label_mtx>0.01] <- ''
png(file.path(dir,'0129-Tissue.LR.sig.png'),res = 600,height = 7,width = 5,units = 'in')
pheatmap(log2(means_df_tumor+1),cluster_rows = F,cluster_cols = F,display_numbers = label_mtx,
         cellwidth = 25,cellheight = 25,scale='row',main = '0129 Tissue')
dev.off()
# ptc
mypvals <- read.table(file.path(dir,"0129_ptc","pvalues.txt"), sep='\t',header = T,check.names = F)
mymeans <- read.table(file.path(dir,"0129_ptc","means.txt"), sep='\t',header = T,check.names = F)
goi <- read.table(file.path(dir,'0129_tumor','GOI.txt'),sep='\t')$V1
coi <- read.table(file.path(dir,'0129_tumor','COI.txt'),sep='\t')$V1

means_df_ptc <- mymeans[mymeans$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(means_df_ptc) <- means_df_ptc$interacting_pair
means_df_ptc <- means_df_ptc[,-1]
pval_df_ptc <- mypvals[mypvals$interacting_pair %in% goi,c('interacting_pair',coi)]
row.names(pval_df_ptc) <- pval_df_ptc$interacting_pair
pval_df_ptc <- pval_df_ptc[,-1]

label_mtx <- as.matrix(pval_df_ptc)
label_mtx[label_mtx<=0.01] <- '\u2217'
label_mtx[label_mtx>0.01] <- ''
png(file.path(dir,'0129-PTC.LR.sig.png'),res = 600,height = 7,width = 5,units = 'in')
pheatmap(log2(means_df_ptc+1),cluster_rows = F,cluster_cols = F,display_numbers = label_mtx,
         cellwidth = 25,cellheight = 25,scale = 'row',main = '0129 PTC')
dev.off()

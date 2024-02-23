library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library("maftools")
library('NMF')
library("MutationalPatterns")
library("ggplot2")
library("ggpubr")
library("tidyverse")
#### LUAD part ######
dir="d:/project_ptc/PTC-lung/mutation/3events_LUAD"
dir_last="d:/project_ptc/PTC-lung/mutation"
files <- lapply(Sys.glob('d:/project_ptc/PTC-lung/mutation/3events_LUAD/*.maf'), read.maf)
tissue_files = files[c(2,4,6,8,10,12,14,16,18,20,22,24)]
ptc_files = files[c(1,3,5,7,9,11,13,15,17,19,21,23)]
merged_maf_luad = merge_mafs(files)
merged_maf_luad_ptc = merge_mafs(ptc_files)
merged_maf_luad_tissue = merge_mafs(tissue_files)

png('d:/project_ptc/PTC-lung/mutation/final_oncoplot.png',res = 600,width = 10,height = 12,units = 'in')
onco_plot = oncoplot(maf = merged_maf_luad,
                     showTitle =F,
                     SampleNamefontSize =1,
                     legendFontSize =1.3,
                     drawColBar =F,
                     drawRowBar=F,
                     draw_titv=F,
                     fill = F,
                     showTumorSampleBarcodes=T,
                     sampleOrder = c('9848-T','9848-C','1113-T','1113-C','1083-T','1083-C','4391-T','4391-C',
                                   'XJZ11-T','XJZ11-C','0289-T','0289-C','XJZ20-T','XJZ20-C','XJZ21-T','XJZ21-C',
                                   'XJZ23-T','XJZ23-C','XJZ25-T','XJZ25-C','XJZ27-T','XJZ27-C','3680-T','3680-C'),
                     top = 40,
                     drawBox =T,
                     fontSize =1,
                     bgCol = '#ffffff')
dev.off()

png(file.path(dir_last,'titv.png'),res = 600,width = 10,height = 12,units = 'in')
merged_maf.titv = titv(maf = merged_maf_luad, plot = F, useSyn = T)
plotTiTv(res = merged_maf.titv,
         showBarcodes =T,
         textSize =1.2,
         sampleOrder = c('9848-T','9848-C','1113-T','1113-C','1083-T','1083-C','4391-T','4391-C',
                         'XJZ11-T','XJZ11-C','0289-T','0289-C','XJZ20-T','XJZ20-C','XJZ21-T','XJZ21-C',
                         'XJZ23-T','XJZ23-C','XJZ25-T','XJZ25-C','XJZ27-T','XJZ27-C','3680-T','3680-C'))
dev.off()


frac = merged_maf.titv$fraction.contribution
png('d:/project_ptc/PTC-lung/mutation/mutationtype.png',res = 600,width = 10,height = 5,units = 'in')
df = data.frame(
  fraction=c(frac$`C>T`,frac$`C>A`,frac$`T>G`,frac$`C>G`,frac$`T>C`,frac$`T>A`),
  sample=rep(c("PTC","Tissue"),72),
  convert=c(rep('C>T',24),rep('C>A',24),rep('T>G',24),rep('C>G',24),rep('T>C',24),rep('T>A',24))
)
p1 = mutate(df,convert = fct_reorder(convert, fraction, .fun='mean' ,.desc = T))%>%
  ggplot(aes(convert, fraction, fill=sample))+geom_boxplot(outlier.shape=NA)+
  theme_minimal()+
  labs(y="%Mutations")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y.left = element_line(size = 1, colour = "grey"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.title.y.left = element_text(size=15))
dev.off()


frac = merged_maf.titv$TiTv.fractions
png('d:/project_ptc/PTC-lung/mutation/titvfraction.png',res = 600,width = 10,height = 5,units = 'in')
df = data.frame(
  fraction=c(frac$Ti,frac$Tv),
  sample=rep(c("PTC","Tissue"),24),
  convert=c(rep('Ti',24),rep('Tv',24))
)
p2=mutate(df,convert = fct_reorder(convert, fraction, .fun='mean' ,.desc = T))%>%
  ggplot(aes(convert, fraction, fill=sample))+geom_boxplot(outlier.shape=NA)+
  theme_minimal()+
  labs(y="%Mutations")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y.left = element_line(size = 1, colour = "grey"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.title.y.left = element_text(size=15))
dev.off()

png('d:/project_ptc/PTC-lung/mutation/titv_final.png',res = 600,width = 15,height = 5,units = 'in')
ggarrange(p1,p2,ncol = 2, nrow = 1)
dev.off()













#### LUSC part ######
dir="d:/project_ptc/PTC-lung/mutation/3events_LUSC"
files <- lapply(Sys.glob('d:/project_ptc/PTC-lung/mutation/3events_LUSC/*.maf'), read.maf)
tissue_files = files[c(2,4,6,8,10,12)]
ptc_files = files[c(1,3,5,7,9,11)]
merged_maf_lusc = merge_mafs(files)
merged_maf_lusc_tissue = merge_mafs(tissue_files)
merged_maf_lusc_ptc = merge_mafs(ptc_files)

png('d:/project_ptc/PTC-lung/mutation/3events_LUSC/final_oncoplot.png',res = 600,width = 10,height = 12,units = 'in')
onco_plot = oncoplot(maf = merged_maf_lusc,
                     showTitle =F,
                     SampleNamefontSize =1,
                     legendFontSize =1.3,
                     drawColBar =F,
                     drawRowBar=F,
                     draw_titv=F,
                     fill = F,
                     showTumorSampleBarcodes=T,
                     sampleOrder = c('2879-T','2879-C','3404-T','3404-C','XJZ10-T','XJZ10-C',
                                     'XJZ29-T','XJZ29-C','XJZ32-T','XJZ32-C',
                                     'XJZ33-T','XJZ33-C'),
                     top = 40,
                     drawBox =T,
                     fontSize =1,
                     bgCol = '#ffffff')
dev.off()

png(file.path(dir,'titv.png'),res = 600,width = 10,height = 12,units = 'in')
merged_maf.titv = titv(maf = merged_maf_lusc, plot = F, useSyn = T)
plotTiTv(res = merged_maf.titv,
         showBarcodes =T,
         textSize =1.2,
         sampleOrder = c('2879-T','2879-C','3404-T','3404-C','XJZ10-T','XJZ10-C',
                         'XJZ29-T','XJZ29-C','XJZ32-T','XJZ32-C',
                         'XJZ33-T','XJZ33-C'))
dev.off()


frac = merged_maf.titv$fraction.contribution
png('d:/project_ptc/PTC-lung/mutation/mutationtype.png',res = 600,width = 10,height = 5,units = 'in')
df = data.frame(
  fraction=c(frac$`C>T`,frac$`C>A`,frac$`T>G`,frac$`C>G`,frac$`T>C`,frac$`T>A`),
  sample=rep(c("PTC","Tissue"),36),
  convert=c(rep('C>T',12),rep('C>A',12),rep('T>G',12),rep('C>G',12),rep('T>C',12),rep('T>A',12))
)
p1=mutate(df,convert = fct_reorder(convert, fraction, .fun='mean' ,.desc = T))%>%
  ggplot(aes(convert, fraction, fill=sample))+geom_boxplot(outlier.shape=NA)+
  theme_minimal()+
  labs(y="%Mutations")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y.left = element_line(size = 1, colour = "grey"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.title.y.left = element_text(size=15))
dev.off()


frac = merged_maf.titv$TiTv.fractions
png('d:/project_ptc/PTC-lung/mutation/titvfraction.png',res = 600,width = 10,height = 5,units = 'in')
df = data.frame(
  fraction=c(frac$Ti,frac$Tv),
  sample=rep(c("PTC","Tissue"),12),
  convert=c(rep('Ti',12),rep('Tv',12))
)
p2=mutate(df,convert = fct_reorder(convert, fraction, .fun='mean' ,.desc = T))%>%
  ggplot(aes(convert, fraction, fill=sample))+geom_boxplot(outlier.shape=NA)+
  theme_minimal()+
  labs(y="%Mutations")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y.left = element_line(size = 1, colour = "grey"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=15),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        axis.title.y.left = element_text(size=15))
dev.off()

png(file.path(dir,'titv_final.png'),res = 600,width = 15,height = 5,units = 'in')
ggarrange(p1,p2,ncol = 2, nrow = 1)
dev.off()

### TCGA cohort ###

png(file.path(dir_last,'tcga.cohort.png'),res = 600,width = 10,height = 10,units = 'in')
laml.mutload = tcgaCompare(maf = c(merged_maf_luad_tissue,merged_maf_luad_ptc,merged_maf_lusc_tissue,merged_maf_lusc_ptc), 
                           cohortName = c('Tissue-LUAD','PTC-LUAD','Tissue-LUSC','PTC-LUSC'),
                           tcga_cohorts =c('LAML','PCPG','THCA','UVM','TGCT','THYM',
                                           'KICH','ACC','LGG','MESO','PRAD','PAAD','BRCA',
                                           'SARC','CHOL','UCS','GBM','KIRC','KIRP','OV',
                                           'UCEC','LIHC','CESC','READ','ESCA','HNSC','DLBC',
                                           'STAD','COAD','LUAD','LUSC','SKCM'),
                           logscale = T, capture_size = 1.1)
dev.off()


#### Pathways ###
png(file.path(dir,'oncogenicPathways.png'),res = 600,width = 15,height = 10,units = 'in')
OncogenicPathways(maf = merged_maf_lusc,
                  fontSize =1.3)
dev.off()
PlotOncogenicPathways(maf = merged_maf_lusc,pathways = c("NOTCH"),
                      showTumorSampleBarcodes = T,fontSize = 1,
                      SampleNamefontSize = 1,
                      sampleOrder = c('2879-T','2879-C','3404-T','3404-C','XJZ10-T','XJZ10-C',
                                      'XJZ29-T','XJZ29-C','XJZ32-T','XJZ32-C',
                                      'XJZ33-T','XJZ33-C'))

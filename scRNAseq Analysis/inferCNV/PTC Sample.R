library(infercnv)
library(Seurat)

# load the Seurat_object
load('combinedCNV_ptc_TNSC.RData')

# code for generating gene_order_file
library(AnnoProbe)
dat<-as.data.frame(GetAssayData(fhs_ptcTN,  slot='counts'))
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
## 这里可以去除性染色体
# 也可以把染色体排序方式改变
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

# visualize the Seurat object
DimPlot(fhs_ptcTN)

# generate infercnv_obj from seurat object
infercnv_obj_full = CreateInfercnvObject(raw_counts_matrix=GetAssayData(fhs_ptcTN,'counts'),
                                         annotations_file=as.matrix(fhs_ptcTN@active.ident),
                                         gene_order_file = geneFile,
                                         delim = '\t',
                                         ref_group_names = c('Fibroblast','Endothelial cell','CAF'))
rm(fhs_ptcTN)
options(scipen = 100)

infercnv_obj_full_run = infercnv::run(infercnv_obj_full,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='PTC_CAF/', 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             leiden_resolution = 0.001,
                             write_expr_matrix = T)

#####
# without CAF as reference cells
infercnv_obj_full_run_noCAF = infercnv::run(infercnv_obj_full,
                                      cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                      out_dir='PTC/', 
                                      cluster_by_groups=TRUE, 
                                      denoise=TRUE,
                                      HMM=TRUE,
                                      leiden_resolution = 0.001)
#####

# plotting
library(tidyverse)
library(vroom)
library(ggpubr)
# get groupings
infercnv_obj = readRDS("PTC_CAF/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data

expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))

colnames(CNV_score)="CNV_score"
CNV_score$Name=rownames(CNV_score)

# load the Seurat obj
load('combinedCNV_ptc_TNSC.RData')
# isolate annotation result
Annotation<-fhs_ptcTN[['annotation']]
Annotation$Name<-rownames(Annotation)
rm(fhs_ptcTN)

# combine
CNV_score<-CNV_score%>%inner_join(Annotation,by='Name')

# plotting
ggboxplot(CNV_score,'annotation','CNV_score',fill = 'annotation',bxp.errorbar = T)+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5,hjust=0.5))

# median score
CNV_score %>%
  group_by(annotation) %>%
  summarise(median_pts = median(CNV_score, na.rm=TRUE))

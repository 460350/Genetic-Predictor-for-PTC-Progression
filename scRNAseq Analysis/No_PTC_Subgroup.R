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

Idents(fhs_ptcTN)<-'overall_ident'

# visualize the Seurat object
DimPlot(fhs_ptcTN)

# generate infercnv_obj from seurat object
infercnv_obj_full = CreateInfercnvObject(raw_counts_matrix=GetAssayData(fhs_ptcTN,'counts'),
                                         annotations_file=as.matrix(fhs_ptcTN@active.ident),
                                         gene_order_file = geneFile,
                                         delim = '\t',
                                         ref_group_names = c('Endothelial cell','Fibroblast'))

rm(fhs_ptcTN)
options(scipen = 100)

# try the optimal analysis_mode subclusters and no pre-determined annotation
infercnv_obj_full_run = infercnv::run(infercnv_obj_full,
                                      cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                      out_dir='PTC_noGrp/', 
                                      cluster_by_groups=F, 
                                      denoise=TRUE,
                                      HMM=TRUE,
                                      leiden_resolution = 0.0001, # default 0.05
                                      analysis_mode = 'subclusters',
                                      tumor_subcluster_pval = 0.01) # default 0.1

# original one
infercnv_obj_full_run = infercnv::run(infercnv_obj_full,
                                      cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                      out_dir='PTC_noSub/', 
                                      cluster_by_groups=TRUE, 
                                      denoise=TRUE,
                                      HMM=TRUE,
                                      leiden_resolution = 0.001)

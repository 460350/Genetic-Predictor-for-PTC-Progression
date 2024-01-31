library(tidyverse)
library(vroom)

# load the normalized expression matrix
exp<-as.data.frame(vroom('thca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt'))
exp<-na.omit(exp)

table(str_sub(colnames(exp[,-c(1,2)]),14,15))
# 01  06 
# 497   1 

# remove the sample with ID ends with 06
m<-str_sub(colnames(exp),15,15)!=6
exp<-exp[,m]
ncol(exp)
table(str_sub(colnames(exp[,-c(1,2)]),14,15))
#  01 
# 497

rownames(exp)<-exp$Hugo_Symbol  # replicate genes: ‘ELMOD1’, ‘FGF13’, ‘NKAIN3’, ‘PALM2AKAP2’, ‘QSOX1(different entrez id)’, ‘SNAP47’, ‘TMEM8B’
# keep the duplicated genes results with higher expression value. Identify them based on order
index<-order(rowMeans(exp[,-c(1,2)]),decreasing = T)
exp_order<-exp[index,]
# remove the duplicated genes results with lower expression value
keep<-!duplicated(exp_order$Hugo_Symbol)
# get the exp without duplicated symbol name
exp<-exp_order[keep,]
rownames(exp)<-exp$Hugo_Symbol

exp<-exp[,-c(1,2)]

# identified the duplicated sample based the sample ID from 1 to 15
l<-!duplicated(str_sub(colnames(exp),1,15));table(l)
# TRUE
# 497


# gene filtration
# only the genes with count > 10 in >50% genes will be kept
k = apply(exp,1, function(x){sum(x>10)>0.5*ncol(exp)});table(k)
# FALSE  TRUE 
# 5956   14555 
exp = exp[k,]
nrow(exp)

# load the clinical metadata
clin<-as.data.frame(vroom('thca_tcga_pan_can_atlas_2018_clinical_data.tsv'))
# clinical info extraction
clin_meta<-clin[,c(
  'Sample ID',
  'Neoplasm Disease Stage American Joint Committee on Cancer Code',
  'Cancer Type Detailed',
  'Disease Free (Months)',
  'Disease Free Status',
  'Months of disease-specific survival',
  'Disease-specific Survival status',
  'Overall Survival (Months)',
  'Overall Survival Status',
  'American Joint Committee on Cancer Metastasis Stage Code',
  'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code',
  'American Joint Committee on Cancer Tumor Stage Code',
  'Person Neoplasm Cancer Status',
  'Progress Free Survival (Months)',
  'Progression Free Status',
  'Primary Lymph Node Presentation Assessment',
  'TMB (nonsynonymous)'
)]
dim(clin_meta)
rownames(clin_meta)<-clin_meta$`Sample ID`
# rename colname
colnames(clin_meta)<-c('ID','Stage','Subtype','DFS_time','DFS_status','DSS_time','DSS_status','OS_time','OS_status',
                       'M_stage','N_stage','T_stage','Tumor_status','PFS_time','PFS_status','Primary_LN_status',
                       'TMB')
# define missing value as NA
clin_meta[clin_meta==""]<-NA

# save the result
save(clin_meta,exp,file = 'WIP_THCA_PanCancer.RData')

table(str_sub(rownames(clin_meta),14,15))
# 01    06 
# 499   1 

# remove the sample with ID ends with 06
n<-str_sub(rownames(clin_meta),15,15)!=6
clin_meta<-clin_meta[n,]
ncol(clin_meta)
table(str_sub(rownames(clin_meta),14,15))
#  01 
# 499

table(clin_meta$PFS_status)
# 0:CENSORED 1:PROGRESSION 
#   448            51

table(clin_meta$DFS_status)
#  0:DiseaseFree 1:Recurred/Progressed 
#     329                    25

table(clin_meta$OS_status)
# 0:LIVING 1:DECEASED 
#  483         16

# modified the stage information
library(stringr)
head(clin_meta$Stage)
table(clin_meta$Stage)
# STAGE I  STAGE II STAGE III  STAGE IV STAGE IVA STAGE IVC 
#  282        51       110         2        46         6 

a = str_extract_all(clin_meta$Stage,"I|V");head(a)
b = sapply(a,paste,collapse = "");head(b)

clin_meta$simp_stage<-b
table(clin_meta$simp_stage)
#  I   II III  IV  NA 
# 282  51 110  54   2 

# remove the samples with NA stage info
clin_meta<-clin_meta[clin_meta$simp_stage!='NA',]

# save the result
save(clin_meta,exp,file = 'WIP_THCA_PanCancer.RData')

# make the clinical metadata based on DFS
clin_meta_DFS<-clin_meta[clin_meta$DFS_status!='NA',]
# Match the clin_data and exp
s = intersect(rownames(clin_meta_DFS),colnames(exp));length(s)
exp_DFS<-exp[,s]
colnames(clin_meta_DFS)
clin_meta_DFS<-clin_meta_DFS[s,]
dim(clin_meta_DFS)

# change the DFS_status
clin_meta_DFS$DFS_status<-ifelse(clin_meta_DFS$DFS_status=='0:DiseaseFree',
                                 0,
                                 1)
table(clin_meta_DFS$DFS_status)
#  0   1 
# 327  25 

# remove the samples without survival time
k1<-clin_meta_DFS$DFS_time>0;table(k1)
# FALSE  TRUE 
#  1     351
clin_meta_DFS<-clin_meta_DFS[k1,]
exp_DFS<-exp_DFS[,k1]

# make the exp matrix to cpm
exp_DFS<-log2(edgeR::cpm(as.matrix(exp_DFS))+1)

# save the result
save(clin_meta_DFS,exp_DFS,file = 'DFS_THCA_PanCancer.RData')

# make the clinical metadata based on PFS
clin_meta_PFS<-clin_meta[clin_meta$PFS_status!='NA',]
# Match the clin_data and exp
s = intersect(rownames(clin_meta_PFS),colnames(exp));length(s)
exp_PFS<-exp[,s]
colnames(clin_meta_DFS)
clin_meta_PFS<-clin_meta_PFS[s,]
dim(clin_meta_PFS)

# change the PFS_status
clin_meta_PFS$PFS_status<-ifelse(clin_meta_PFS$PFS_status=='0:CENSORED',
                        0,
                        1)
table(clin_meta_PFS$PFS_status)
#  0   1 
# 446  51 

# remove the samples without survival time
k2<-clin_meta_PFS$PFS_time>0;table(k2)
# FALSE  TRUE 
#  1     496
clin_meta_PFS<-clin_meta_PFS[k2,]
exp_PFS<-exp_PFS[,k2]

# make the exp matrix to cpm
exp_PFS<-log2(edgeR::cpm(as.matrix(exp_PFS))+1)

# save the result
save(clin_meta_PFS,exp_PFS,file = 'PFS_THCA_PanCancer.RData')

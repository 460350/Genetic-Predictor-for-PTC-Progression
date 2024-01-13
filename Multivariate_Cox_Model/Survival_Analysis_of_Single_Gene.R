# DATA IMPORT
rm(list=ls())
project = "TCGA-KIRC"
load(paste0("data/",project,"_gdc.Rdata"))
library(stringr)

# EXPRESSION DATA PREPARATION
## remove normal sample in expression matrix
table(Group)

exprSet<-exp[,Group=='tumor']
ncol(exp)
ncol(exprSet)

## gene filtration
# only the genes with count > 10 in >50% genes will be kept
k = apply(exprSet,1, function(x){sum(x>10)>0.5*ncol(exprSet)});table(k)
# FALSE  TRUE 
# 11265 19672 
exprSet = exprSet[k,]
nrow(exprSet)

# cpm calculation
cpm_exprSet<-log2(edgeR::cpm(exprSet)+1)

# remove the duplicate sample
# reorder the expression matrix based on the sample ID
cpm_exprSet<-cpm_exprSet[,sort(colnames(cpm_exprSet))]
# identified the duplicated sample based the sample ID from 1 to 12
l<-!duplicated(str_sub(colnames(cpm_exprSet),1,12));table(l)
# FALSE  TRUE 
#  8     503 
cpm_exprSet<-cpm_exprSet[,l]
colnames(cpm_exprSet)=str_sub(colnames(cpm_exprSet),1,12)
ncol(cpm_exprSet)

# CLINICAL INFORMATION PREPARATION
# extract wanted info
clin_meta<-clinical[,c(
  'bcr_patient_barcode',
  'vital_status',
  'days_to_death',
  'days_to_last_followup',
  'race_list',
  'age_at_initial_pathologic_diagnosis',
  'gender',
  'stage_event',
  "braf_gene_genotyping_outcome_lab_results_text",
  "ras_family_gene_genotyping_outcome_lab_results_text",
  "ret_ptc_rearrangement_genotyping_outcome_lab_results_text",
  "new_tumor_events"
)]
dim(clin_meta)
rownames(clin_meta)<-clin_meta$bcr_patient_barcode
# rename colname
colnames(clin_meta)<-c('ID','event','death','last_followup','race','age','gender','stage',
                       'BRAF_mut','RAS_mut',"RET_rearrange","recurrence")
# define missing value as NA
clin_meta[clin_meta==""]<-NA

# calculation the following-up time
table(clin_meta$event)
# Alive  Dead 
#  493    14 
clin_meta$time<-ifelse(clin_meta$event=='Alive',
                       clin_meta$last_followup,
                       clin_meta$death)
clin_meta$time<-as.numeric(clin_meta$time)/30
clin_meta$age<-as.numeric(clin_meta$age)

# define event, Alive = 0, Dead = 1
clin_meta$event<-ifelse(clin_meta$event=='Alive',
                        0,
                        1)
table(clin_meta$event)
#  0   1 
# 493  14 

# remove the samples without survival time
k1<-clin_meta$time>0.1;table(k1)
# FALSE  TRUE 
#  13    493 
# remove the samples without survival time or dead/alive event
k2<-!(is.na(clin_meta$time)|is.na(clin_meta$event));table(k2)
# FALSE  TRUE 
#   1    506 

clin_meta<-clin_meta[k1&k2,]
nrow(clin_meta)
# 493

# modified the stage information
library(stringr)
head(clin_meta$stage)

a = str_extract_all(clin_meta$stage,"I|V");head(a)
b = sapply(a,paste,collapse = "");head(b)

clin_meta$simp_stage<-b
table(clin_meta$simp_stage)

#    I   II III  IV 
# 2 275  52 110  54 

# set the missing value as NA
clin_meta[clin_meta==""]<-as.character(NA)
table(clin_meta$simp_stage,useNA = "always")
#   I    II  III   IV  <NA> 
#  275   52  110   54    2 

# Match the clin_data and cpm_exprSet
s = intersect(rownames(clin_meta),colnames(cpm_exprSet));length(s)
cpm_exprSet<-cpm_exprSet[,s]
colnames(clin_meta)
clin_meta<-clin_meta[s,c(1,2,13,4:7,14,9:12)]
dim(clin_meta)

# SURVIVAL ANALYSIS OF SINGLE GENE
# analysis the impact of single gene expression on the survival
genelist<-c('FN1','LGALS3','CAV1','APOE')
library(survival)
library(survminer)
splots = list()

for (i in 1:length(genelist)) {
  x = clin_meta
  g = genelist[i]
  x$gene = ifelse(cpm_exprSet[g,]>median(cpm_exprSet[g,]),'high','low')
  sfit1 = survfit(Surv(time,event)~gene, data = x)
  splots[[i]] = ggsurvplot(sfit1, pval = T, palette = c("red","blue"),
                           data = x, legend = c(0.8,0.8),
                           title = genelist[i],risk.table = T)
}
splots

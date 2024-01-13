library(survival)
library(survminer)
library(dplyr)
library(glmnet)
library(tidyverse)
library(GGally)
library(rms)
library(survivalROC)
library(plotROC)
library(stringr)
library(vroom)

# identify the significant genes and give them order based on padj
resOrder<-res[order(res$padj),]
resSig<-subset(resOrder,padj<0.01)
resSigAll<-subset(resSig,(log2FoldChange<(-2)|log2FoldChange>2))
# normalize the data with vst approach
rld<-vst(dds,blind = T)
expr_norm<-assay(rld)
DESeq_norm_vst_for_survival<-expr_norm[resSigAll@rownames,]

# prepare the survival table
survival_cancer<-clin_meta%>%select(time,event,simp_stage)
# in package 'survival', 0 represents nothing happens, 1 represents diseased/dead
colnames(survival_cancer)<-c("overall_survival","censoring_status","tumor_stage")
survival_cancer$overall_survival<-as.numeric(survival_cancer$overall_survival)
survival_cancer$censoring_status<-as.numeric(survival_cancer$censoring_status)

# remove the NA column
survival_cancer<-survival_cancer[!is.na(survival_cancer),]
rownames(survival_cancer)
# keep the rows with real data
survival_cancer<-survival_cancer[1:489,]
survival_cancer$TCGA_IDs<-rownames(survival_cancer)
rownames(survival_cancer)

# check the rownames again with m.Targeted.list
survival_cancer<-survival_cancer%>%
  dplyr::filter(TCGA_IDs %in% rownames(clin_meta))
# remove the rows with overall_survival=0
survival_cancer<-survival_cancer%>%
  dplyr::filter(overall_survival != 0)

# remove the normal sample among DESeq_norm_vst_for_survival and change the name from 15-digit to 12-digit
# remove the normal sample
table(str_sub(colnames(DESeq_norm_vst_for_survival),14,15))
t<-str_sub(colnames(DESeq_norm_vst_for_survival),14,14)!=1
DESeq_norm_vst_for_survival<-DESeq_norm_vst_for_survival[,t]
table(str_sub(colnames(DESeq_norm_vst_for_survival),14,15))
# change the name from 15-digit to 12 digit
colnames(DESeq_norm_vst_for_survival)=str_sub(colnames(DESeq_norm_vst_for_survival),1,12)

# add the expression data to survival_cancer
survival_cancer<-cbind(survival_cancer,t(DESeq_norm_vst_for_survival)[survival_cancer$TCGA_IDs,])
# change the gene names containing '-'
colnames(survival_cancer)<-gsub(colnames(survival_cancer),pattern='-',replacement='_')

# save the survival data and expression data with cancer-only data
save(survival_cancer,clin_meta,file = 'TCGA_THCA_phenotype&expr.RData')

# univariate cox regression
# filter potential useful sig genes using univariate cox regression
uni_cox_in_bulk<-function(gene_list,survival_info_df){
  library(survival)
  gene_list<-gsub(gene_list,pattern = '-',replacement = '_')
  uni_cox<-function(single_gene){
    formula<-as.formula(paste0('Surv(overall_survival,censoring_status)~',single_gene))
    surv_uni_cox<-summary(coxph(formula,data=survival_cancer))
    ph_hypothesis_p<-cox.zph(coxph(formula,data = survival_cancer))$table[1,3]
    if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){ # get the pvalue
      single_cox_report<-data.frame("uni_cox_sig_genes"=single_gene,
                                    "beta"=surv_uni_cox$coefficients[,1],
                                    "Hazard_Ratio"=exp(surv_uni_cox$coefficients[,1]),
                                    "z_pvalue"=surv_uni_cox$coefficients[,5],
                                    "Wald_pvalue"=as.numeric(surv_uni_cox$waldtest[3]),
                                    "Likelihood_pvalue"=as.numeric(surv_uni_cox$logtest[3]))
      single_cox_report
    }
  }
  uni_cox_list<-lapply(gene_list,uni_cox)
  do.call(rbind,uni_cox_list)
}
uni_cox_df<-uni_cox_in_bulk(gene_list = resSigAll@rownames,survival_info_df = survival_cancer)
uni_cox_sig_genes<-uni_cox_df$uni_cox_sig_genes
uni_cox_sig_genes
# 140 genes identified
# ZCCHC12 is one of them 

# use Lasso regression
x<-as.matrix(survival_cancer[,gsub(resSigAll@rownames,
                                   pattern = '-',replacement='_')])
y<-survival_cancer[,c('overall_survival','censoring_status')]
names(y)<-c('time','status')
y$time<-as.double(y$time)
y$status<-as.double(y$status)
y<-as.matrix(survival::Surv(y$time,y$status))
lasso_fit<-cv.glmnet(x,y,family='cox',type.measure = 'deviance')
coefficient<-coef(lasso_fit,s=lasso_fit$lambda.min)
Active.Index<-which(as.numeric(coefficient) != 0)
active.coefficients<-as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox<-rownames(coefficient)[Active.Index]
sig_gene_multi_cox
# 20 genes identified
# EEF1A2 is one of them

# use the candidates identifed above to construct multivariate model
# perform the multi-variates cox regression using qualified genes
formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(sig_gene_multi_cox,sep='',collapse = '+')))
multi_variate_cox<-coxph(formula_for_multivariate,data = survival_cancer)
# check if variances are supported by PH hypothesis
ph_hypo_multi<-cox.zph(multi_variate_cox)
# The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table<-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
# remove variances not supported by ph hypothesis and perform the 2nd regression
formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],
                                                  sep='',collapse = '+')))
multi_variate_cox_2<-coxph(formula_for_multivariate,data = survival_cancer)
ph_hypo_multi_2<-cox.zph(multi_variate_cox_2)
ph_hypo_table_2<-ph_hypo_multi_2$table[-nrow(ph_hypo_multi_2$table),]

# check the co-linearity between samples
correlation<-cor(survival_cancer[,rownames(ph_hypo_table_2)[ph_hypo_table_2[,3]>0.05]],
                 method = 'pearson')
library(GGally)
ggpairs(survival_cancer[,rownames(ph_hypo_table_2)[ph_hypo_table_2[,3]>0.05]],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())
library(rms)
vif<-rms::vif(multi_variate_cox_2)
sqrt(vif)<2
# keep the true result
# 'RIPPLY3', 'HPSE2', 'DUXAP8', 'CNGA3', 'NLGN1', 'GRIN1', 'ZSCAN4', 'APOD'
ggpairs(survival_cancer[,c('RIPPLY3','HPSE2','DUXAP8','CNGA3','NLGN1','GRIN1','ZSCAN4','APOD')],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())

formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(c('RIPPLY3','HPSE2','DUXAP8','CNGA3','NLGN1','GRIN1','ZSCAN4','APOD'),
                                                  sep='',collapse = '+')))
multi_variate_cox_3<-coxph(formula_for_multivariate,data = survival_cancer)
ph_hypo_multi_3<-cox.zph(multi_variate_cox_3)
ph_hypo_table_3<-ph_hypo_multi_3$table[-nrow(ph_hypo_multi_3$table),]


# forest plot
# HR<1 for protection, HR>1 for risk
ggforest(model = multi_variate_cox_3,
         data = survival_cancer,
         main = "Hazard ratios of candidate genes",
         fontsize = 1)

# evaluate the cox model
## concordance index
C_index<-multi_variate_cox_3$concordance['concordance']
if(C_index>=0.9){
  print('High accuracy')
}else{
  if(C_index<0.9&C_index>=0.7){
    print('Medium accuracy')
  }else{
    print('Low accuracy')
  }
}
# "High accuracy"

## time-dependent ROC curve
### usually AUC>=0.7 is good
# calculate the risk score of each sample
riskscore<-function(survival_cancer_df, candidate_genes_for_cox, cox_report){
  library(dplyr)
  risk_score_table<-survival_cancer_df[,candidate_genes_for_cox]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene<-risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table<-cbind(risk_score_table,'total_risk_score'=exp(rowSums(risk_score_table)))%>%
    cbind(survival_cancer_df[,c('TCGA_IDs','overall_survival','censoring_status')])
  risk_score_table<-risk_score_table[,c('TCGA_IDs','overall_survival','censoring_status',candidate_genes_for_cox,'total_risk_score')]
  risk_score_table
}
candidate_genes_for_cox3<-c(rownames(ph_hypo_table_3)[ph_hypo_table_3[,3]>0.05])
risk_score_table_multi_cox3<-riskscore(survival_cancer,
                                       candidate_genes_for_cox3,
                                       multi_variate_cox_3)

# make ROC curve
multi_ROC<-function(time_vector,risk_score_table){
  library(survivalROC)
  single_ROC<-function(single_time){
    for_ROC<-survivalROC(Stime = risk_score_table$overall_survival,
                         status = risk_score_table$censoring_status,
                         marker = risk_score_table$total_risk_score,
                         predict.time = single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP,'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values,'Time_point'=rep(single_time,length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC,length(for_ROC$TP)))
  }
  multi_ROC_list<-lapply(time_vector,single_ROC)
  do.call(rbind,multi_ROC_list)
}
# evaluate 11 AUCs between 3-5 years
for_multi_ROC<-multi_ROC(time_vector = c(365*seq(3,5,0.2)),
                         risk_score_table = risk_score_table_multi_cox3)
# multiple lines are overlapped

# KM-curve
# use the time with max AUC to do KM analysis
AUC_max<-max(for_multi_ROC$AUC)
# select the last time point with max AUC indicating longer survival
AUC_max_time<-for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
# remove duplicate
AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
AUC_max_time<-AUC_max_time[length(AUC_max_time)]
for_multi_ROC$Time_point<-as.factor(for_multi_ROC$Time_point)
# find the optimal cutoff value within the ROC curve of the optimal time point
optimal_time_ROC_df<-for_multi_ROC[which(for_multi_ROC$Time_point==AUC_max_time),]
cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive-optimal_time_ROC_df$False_positive)]
high_low<-(risk_score_table_multi_cox3$total_risk_score>cut.off)
high_low[high_low==T]<-'high'
high_low[high_low==F]<-'low'
risk_score_table_multi_cox3<-cbind(risk_score_table_multi_cox3,high_low)

# KM_plot generation
library(survminer)
# first edit the status of patients with OS > AUC max time. (censoring status=0 (alive), OS=365*5 days)
risk_score_table_multi_cox3$censoring_status[which(risk_score_table_multi_cox3$overall_survival>AUC_max_time)]<-0
risk_score_table_multi_cox3$overall_survival[which(risk_score_table_multi_cox3$overall_survival>AUC_max_time)]<-AUC_max_time
fit_km<-survfit(Surv(overall_survival,censoring_status)~high_low,data = risk_score_table_multi_cox3)
ggsurvplot(fit_km,conf.int = F,pval = T,legend.title='total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),surv.median.line = 'hv')

# visualization of the ROC curves of multiple time point
pROC<-ggplot(for_multi_ROC,
             aes(x=False_positive,y=True_positive,label=Cut_values,color=Time_point))+
  geom_roc(labels = F,stat = 'identity',n.cuts = 0)+
  geom_abline(slope = 1,intercept = 0,color='red',linetype=2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill='white'),
        panel.grid = element_blank())+
  annotate("text",x=0.75,y=0.15,
           label=paste("AUC max = ",round(AUC_max,2),'\n', 'AUC max time =', 
                       AUC_max_time, ' days', sep = ''))
pROC
# AUCmax=0.59
# THCA has good prognostic, the KM plot may not be good enough to differentiate the result?


# construct model from uni_cox_sig_genes
formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(uni_cox_sig_genes,sep='',collapse = '+')))
multi_variate_cox_4<-coxph(formula_for_multivariate,data = survival_cancer)
# it fails to converge
# check if variances are supported by PH hypothesis
ph_hypo_multi_4<-cox.zph(multi_variate_cox_4)
# system is computationally singular: reciprocal condition number = 2.10145e-16

# save the essential compartment for analysis using other genelists
# it inlcude expr_norm of the TCGA_THCA, survival_cancer data.frame without the expression matrix
# the process here have already removed the normal sample in the expression matrix,
# aligned the sample names between survival data and expression matrix, 
# and removed the NA value in the survival_cancer
surv_ca<-survival_cancer[,1:4]
save(expr_norm,surv_ca,file = paste0("data/Basic_Data_for_TCGA_THCA_Survival_Analysis.Rdata"))

# try the DEGs from Scissor analysis
Scissor.de.gene<-as.data.frame(vroom('Scissor_stage_v2_ptc_de_genes.csv'))
Scissor.de.gene<-Scissor.de.gene[,-1]
rownames(Scissor.de.gene)<-Scissor.de.gene$SYMBOL

# remove the genes that are not included in the data.frame survival cancer
Scissor.de.gene<-Scissor.de.gene%>%
  dplyr::filter(SYMBOL %in% colnames(survival_cancer[,-c(1:4)]))

# univariate cox regression
# filter potential useful sig genes using univariate cox regression
uni_cox_df<-uni_cox_in_bulk(gene_list = Scissor.de.gene$SYMBOL,survival_info_df = survival_cancer)
uni_cox_sig_genes<-uni_cox_df$uni_cox_sig_genes
uni_cox_sig_genes
# 6 genes identified

# use Lasso regression
x<-as.matrix(survival_cancer[,gsub(Scissor.de.gene$SYMBOL,
                                   pattern = '-',replacement='_')])
y<-survival_cancer[,c('overall_survival','censoring_status')]
names(y)<-c('time','status')
y$time<-as.double(y$time)
y$status<-as.double(y$status)
y<-as.matrix(survival::Surv(y$time,y$status))
lasso_fit<-cv.glmnet(x,y,family='cox',type.measure = 'deviance')
coefficient<-coef(lasso_fit,s=lasso_fit$lambda.min)
Active.Index<-which(as.numeric(coefficient) != 0)
active.coefficients<-as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox<-rownames(coefficient)[Active.Index]
sig_gene_multi_cox
# 20 genes identified

# use the candidates identifed above to construct multivariate model
# perform the multi-variates cox regression using qualified genes
formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(sig_gene_multi_cox,sep='',collapse = '+')))
multi_variate_cox<-coxph(formula_for_multivariate,data = survival_cancer)
# check if variances are supported by PH hypothesis
ph_hypo_multi<-cox.zph(multi_variate_cox)
# The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table<-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
# remove variances not supported by ph hypothesis and perform the 2nd regression
formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],
                                                  sep='',collapse = '+')))
multi_variate_cox_2<-coxph(formula_for_multivariate,data = survival_cancer)
ph_hypo_multi_2<-cox.zph(multi_variate_cox_2)
ph_hypo_table_2<-ph_hypo_multi_2$table[-nrow(ph_hypo_multi_2$table),]

# check the co-linearity between samples
correlation<-cor(survival_cancer[,rownames(ph_hypo_table_2)[ph_hypo_table_2[,3]>0.05]],
                 method = 'pearson')
library(GGally)
ggpairs(survival_cancer[,rownames(ph_hypo_table_2)[ph_hypo_table_2[,3]>0.05]],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())
library(rms)
vif<-rms::vif(multi_variate_cox_2)
vif<5

# keep the true result
# 'APOE', 'NMU', 'CRABP2', 'PNP', 'PASK', 'CBLN1'
ggpairs(survival_cancer[,c('APOE', 'NMU', 'CRABP2', 'PNP', 'PASK', 'CBLN1')],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())

formula_for_multivariate<-as.formula(paste0('Surv(overall_survival,censoring_status)~',
                                            paste(c('APOE', 'NMU', 'CRABP2', 'PNP', 'PASK', 'CBLN1'),
                                                  sep='',collapse = '+')))
multi_variate_cox_3<-coxph(formula_for_multivariate,data = survival_cancer)
ph_hypo_multi_3<-cox.zph(multi_variate_cox_3)
ph_hypo_table_3<-ph_hypo_multi_3$table[-nrow(ph_hypo_multi_3$table),]

# forest plot
# HR<1 for protection, HR>1 for risk
ggforest(model = multi_variate_cox_2,
         data = survival_cancer,
         main = "Hazard ratios of candidate genes",
         fontsize = 1)

# evaluate the cox model
## concordance index
C_index<-multi_variate_cox_2$concordance['concordance']
if(C_index>=0.9){
  print('High accuracy')
}else{
  if(C_index<0.9&C_index>=0.7){
    print('Medium accuracy')
  }else{
    print('Low accuracy')
  }
}
# "High accuracy"

## time-dependent ROC curve
### usually AUC>=0.7 is good
# calculate the risk score of each sample
candidate_genes_for_cox2<-c(rownames(ph_hypo_table_2)[ph_hypo_table_2[,3]>0.05])
risk_score_table_multi_cox2<-riskscore(survival_cancer,
                                       candidate_genes_for_cox2,
                                       multi_variate_cox_2)

# make ROC curve
# evaluate 11 AUCs between 3-5 years
for_multi_ROC<-multi_ROC(time_vector = c(365*seq(3,5,0.2)),
                         risk_score_table = risk_score_table_multi_cox2)
# multiple lines are overlapped

# KM-curve
# use the time with max AUC to do KM analysis
AUC_max<-max(for_multi_ROC$AUC)
# select the last time point with max AUC indicating longer survival
AUC_max_time<-for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
# remove duplicate
AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
AUC_max_time<-AUC_max_time[length(AUC_max_time)]
for_multi_ROC$Time_point<-as.factor(for_multi_ROC$Time_point)
# find the optimal cutoff value within the ROC curve of the optimal time point
optimal_time_ROC_df<-for_multi_ROC[which(for_multi_ROC$Time_point==AUC_max_time),]
cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive-optimal_time_ROC_df$False_positive)]
high_low<-(risk_score_table_multi_cox2$total_risk_score>cut.off)
high_low[high_low==T]<-'high'
high_low[high_low==F]<-'low'
risk_score_table_multi_cox2<-cbind(risk_score_table_multi_cox2,high_low)

# KM_plot generation
library(survminer)
# first edit the status of patients with OS > AUC max time. (censoring status=0 (alive), OS=365*5 days)
risk_score_table_multi_cox2$censoring_status[which(risk_score_table_multi_cox2$overall_survival>AUC_max_time)]<-0
risk_score_table_multi_cox2$overall_survival[which(risk_score_table_multi_cox2$overall_survival>AUC_max_time)]<-AUC_max_time
fit_km<-survfit(Surv(overall_survival,censoring_status)~high_low,data = risk_score_table_multi_cox2)
ggsurvplot(fit_km,conf.int = F,pval = T,legend.title='total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),surv.median.line = 'hv')

# visualization of the ROC curves of multiple time point
pROC<-ggplot(for_multi_ROC,
             aes(x=False_positive,y=True_positive,label=Cut_values,color=Time_point))+
  geom_roc(labels = F,stat = 'identity',n.cuts = 0)+
  geom_abline(slope = 1,intercept = 0,color='red',linetype=2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill='white'),
        panel.grid = element_blank())+
  annotate("text",x=0.75,y=0.15,
           label=paste("AUC max = ",round(AUC_max,2),'\n', 'AUC max time =', 
                       AUC_max_time, ' days', sep = ''))
pROC
# AUCmax=0.45
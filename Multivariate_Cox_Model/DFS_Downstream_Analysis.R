library(tidyverse)
library(vroom)
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

# load the DFS data
load('DFS_THCA_PanCancer.RData')

# load the common genes from Scissor DFS&stage V2 analysis
GeneList<-as.data.frame(vroom('Combination_of_DFS_Stage_V2.csv',delim = "/t"))
# make it a list and remove the duplicate genes
C_list<-GeneList$Gene

k3<-!duplicated(C_list)
C_list<-C_list[k3]

####COX_KM_ANALYSIS####
# prepare the survival table
survival_cancer<-clin_meta_DFS%>%select(DFS_time,DFS_status,simp_stage)
# in package 'survival', 0 represents nothing happens, 1 represents diseased/dead
colnames(survival_cancer)<-c("DFS_time","DFS_status","tumor_stage")
survival_cancer$DFS_time<-as.numeric(survival_cancer$DFS_time)
survival_cancer$DFS_status<-as.numeric(survival_cancer$DFS_status)

# change the simp_stage as number
survival_cancer$tumor_stage[survival_cancer$tumor_stage == 'I'] <- 1
survival_cancer$tumor_stage[survival_cancer$tumor_stage == 'II'] <- 2
survival_cancer$tumor_stage[survival_cancer$tumor_stage == 'III'] <- 3
survival_cancer$tumor_stage[survival_cancer$tumor_stage == 'IV'] <- 4
survival_cancer$tumor_stage<-as.numeric(survival_cancer$tumor_stage)
str(survival_cancer)

table(is.na(survival_cancer))

# add the expression data to survival_cancer
survival_cancer$TCGA_IDs<-rownames(survival_cancer)
rownames(survival_cancer)
survival_cancer<-cbind(survival_cancer,t(exp_DFS)[survival_cancer$TCGA_IDs,])
# change the colnames of survival_cancer
colnames(survival_cancer)<-gsub(colnames(survival_cancer),pattern = '-',replacement = '_')

# evaluate the relationship between simp_stage and DFS status
stage_DFSstatus<-clin_meta_DFS[,c(5,18)]
# change I, II, III into 1, 2, 3
stage_DFSstatus$simp_stage[stage_DFSstatus$simp_stage=='I']<-1
stage_DFSstatus$simp_stage[stage_DFSstatus$simp_stage=='II']<-2
stage_DFSstatus$simp_stage[stage_DFSstatus$simp_stage=='III']<-3
stage_DFSstatus$simp_stage<-as.numeric(stage_DFSstatus$simp_stage)
table(stage_DFSstatus)

# evaluate the correlation
shapiro.test(stage_DFSstatus$DFS_status)
# p-value < 2.2e-16
shapiro.test(stage_DFSstatus$simp_stage)
# p-value < 2.2e-16
# data are not normally distributed, using kendall instead of pearson
cor_stageDFS<-cor(stage_DFSstatus$DFS_status,stage_DFSstatus$simp_stage,method = 'kendall')
cor_stageDFS
# 0.070965
cor.test(stage_DFSstatus$DFS_status,stage_DFSstatus$simp_stage,method = 'kendall')

# save the result
save(clin_meta_DFS,exp_DFS,survival_cancer,file = 'DFS_THCA_PanCancer.RData')

# univariate cox regression
# filter potential useful sig genes using univariate cox regression
uni_cox_in_bulk<-function(gene_list,survival_info_df){
  library(survival)
  gene_list<-gsub(gene_list,pattern = '-',replacement = '_')
  uni_cox<-function(single_gene){
    formula<-as.formula(paste0('Surv(DFS_time,DFS_status)~',single_gene))
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
# add tumor_stage to the C_list
C_list<-append('tumor_stage',C_list)
# remove the genes that not among the colnames of survival cancer
C_list<-C_list[C_list%in%colnames(survival_cancer)]

uni_cox_df<-uni_cox_in_bulk(gene_list = C_list,survival_info_df = survival_cancer)
uni_cox_sig_genes<-uni_cox_df$uni_cox_sig_genes
uni_cox_sig_genes
# 2 genes are identified, HMGB2 & STMN1

# use Lasso regression
x<-as.matrix(survival_cancer[,gsub(C_list,
                                   pattern = '-',replacement='_')])
y<-survival_cancer[,c('DFS_time','DFS_status')]
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
# 0 gene is identified

# try all the genes instead
A_list<-rownames(exp_DFS)
# add tumor_stage to the A_list
A_list<-append('tumor_stage',A_list)
uni_cox_df<-uni_cox_in_bulk(gene_list = A_list,survival_info_df = survival_cancer)
uni_cox_sig_genes<-uni_cox_df$uni_cox_sig_genes
uni_cox_sig_genes
# 1321 genes have been identified

# use Lasso regression
x<-as.matrix(survival_cancer[,gsub(A_list,
                                   pattern = '-',replacement='_')])
y<-survival_cancer[,c('DFS_time','DFS_status')]
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
# 39 gene are identified

# use the candidates identifed above to construct multivariate model
# perform the multi-variates cox regression using qualified genes
formula_for_multivariate<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                            paste(sig_gene_multi_cox,sep='',collapse = '+')))
multi_variate_cox<-coxph(formula_for_multivariate,data = survival_cancer)
# check if variances are supported by PH hypothesis
ph_hypo_multi<-cox.zph(multi_variate_cox)
# The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table<-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
# remove variances not supported by ph hypothesis and perform the 2nd regression
formula_for_multivariate<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
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
sqrt(vif)<5
rms_T<-names(ph_hypo_table_2[,3][sqrt(vif)<5])

ggpairs(survival_cancer[,rms_T],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())

formula_for_multivariate<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                            paste(rms_T,
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
    cbind(survival_cancer_df[,c('TCGA_IDs','DFS_time','DFS_status')])
  risk_score_table<-risk_score_table[,c('TCGA_IDs','DFS_time','DFS_status',candidate_genes_for_cox,'total_risk_score')]
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
    for_ROC<-survivalROC(Stime = risk_score_table$DFS_time,
                         status = risk_score_table$DFS_status,
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
risk_score_table_multi_cox3$DFS_status[which(risk_score_table_multi_cox3$DFS_time>AUC_max_time)]<-0
risk_score_table_multi_cox3$DFS_time[which(risk_score_table_multi_cox3$DFS_time>AUC_max_time)]<-AUC_max_time
fit_km<-survfit(Surv(DFS_time,DFS_status)~high_low,data = risk_score_table_multi_cox3)
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
# AUCmax=0.54


####directly use C_list to construct multi variate cox model####

# use the candidates identifed above to construct multivariate model
# perform the multi-variates cox regression using qualified genes
formula_for_multivariate_c<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                            paste(gsub(C_list,
                                                       pattern = '-',replacement='_'),sep='',collapse = '+')))
multi_variate_cox_c<-coxph(formula_for_multivariate_c,data = survival_cancer)
# check if variances are supported by PH hypothesis
ph_hypo_multi_c<-cox.zph(multi_variate_cox_c)
# The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table_c<-ph_hypo_multi_c$table[-nrow(ph_hypo_multi_c$table),]
# remove variances not supported by ph hypothesis and perform the 2nd regression
formula_for_multivariate_c<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                            paste(rownames(ph_hypo_table_c)[ph_hypo_table_c[,3]>0.05],
                                                  sep='',collapse = '+')))
multi_variate_cox_c_2<-coxph(formula_for_multivariate_c,data = survival_cancer)
ph_hypo_multi_c_2<-cox.zph(multi_variate_cox_c_2)
ph_hypo_table_c_2<-ph_hypo_multi_c_2$table[-nrow(ph_hypo_multi_c_2$table),]

# check the co-linearity between samples
correlation_c<-cor(survival_cancer[,rownames(ph_hypo_table_c_2)[ph_hypo_table_c_2[,3]>0.05]],
                 method = 'pearson')
library(GGally)
ggpairs(survival_cancer[,rownames(ph_hypo_table_c_2)[ph_hypo_table_c_2[,3]>0.05]],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())
library(rms)
vif_c<-rms::vif(multi_variate_cox_c_3)
sqrt(vif_c)<3

# some variable the colinearirity is significant, reduce it with lasso
# define response variable
Y<-as.matrix(survival::Surv(survival_cancer$DFS_time,survival_cancer$DFS_status))

# define a matrix of predictor variables
X<-as.matrix(survival_cancer[,gsub(rownames(ph_hypo_table_c_2),
                                   pattern = '-',replacement='_')])

fit<-glmnet(X,Y,family = 'cox',alpha=1)
plot(fit,label = T)
plot(fit,xvar="lambda",label=T)

# k-fold cross-validation to fine the optimal lambda value
fitcv <- cv.glmnet(X,Y,family="cox", alpha=1,nfolds=10)
# visualization of the optimal lambda value that minimizes the result
plot(fitcv)
coef(fitcv, s="lambda.min")
best_lambda<-fitcv$lambda.min
best_lambda

# alternative starting fromt the cross validation
# default type.measure="deviance", here try C-index.Both are good for cox model.
fitcv_test <- cv.glmnet(X,Y,family="cox", alpha=1,nfolds=10,type.measure = 'C')
plot(fitcv_test)
coef(fitcv_test, s="lambda.min")
best_lambda_test<-fitcv_test$lambda.min
best_lambda_test
# result in  14 or 10 variables, choose the result of deviance instead

# test the best model
best_model<-glmnet(X,Y,family='cox',alpha=1,lambda = best_lambda)
coef(best_model)
print(best_model)

# demonstrate the cox model
coxm<-cph(Surv(DFS_time,DFS_status==1)~MAPK13+HMGB2+STMN1,x=T,y=T,data=survival_cancer,surv=T)
cox.zph(coxm)
#        chisq df   p
# MAPK13 0.926  1 0.34
# HMGB2  2.215  1 0.14
# STMN1  0.686  1 0.41
# GLOBAL 3.845  3 0.28
coxm
plot(cox.zph(coxm))
ggcoxzph(cox.zph(coxm))
# ideally the cox plot should be a line when k=0, suggesting that y won't change as time changes
ggcoxdiagnostics(coxm, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# Specifying the argument type = “dfbeta”, plots the estimated changes in the regression coefficients upon deleting each observation in turn
ggcoxdiagnostics(coxm, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# 'The pattern looks fairly symmetric around 0' is the best result.
ggcoxfunctional(formula_for_multivariate_c,data = survival_cancer)
# MAPK13, HMGB2, STMN1 are identified

ggpairs(survival_cancer[,c('MAPK13','HMGB2','STMN1')],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())

formula_for_multivariate_c<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                            paste(c('MAPK13','HMGB2','STMN1'),
                                                  sep='',collapse = '+')))
multi_variate_cox_c_3<-coxph(formula_for_multivariate_c,data = survival_cancer)
ph_hypo_multi_c_3<-cox.zph(multi_variate_cox_c_3)
ph_hypo_table_c_3<-ph_hypo_multi_c_3$table[-nrow(ph_hypo_multi_c_3$table),]
plot(cox.zph(multi_variate_cox_c_3))
ggcoxzph(cox.zph(multi_variate_cox_c_3))

# forest plot
# HR<1 for protection, HR>1 for risk
ggforest(model = multi_variate_cox_c_3,
         data = survival_cancer,
         main = "Hazard ratios of candidate genes",
         fontsize = 1)

# evaluate the cox model
## concordance index
C_index<-multi_variate_cox_c_3$concordance['concordance']
if(C_index>=0.9){
  print('High accuracy')
}else{
  if(C_index<0.9&C_index>=0.7){
    print('Medium accuracy')
  }else{
    print('Low accuracy')
  }
}
C_index
# "Medium accuracy" 0.74

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
    cbind(survival_cancer_df[,c('TCGA_IDs','DFS_time','DFS_status')])
  risk_score_table<-risk_score_table[,c('TCGA_IDs','DFS_time','DFS_status',candidate_genes_for_cox,'total_risk_score')]
  risk_score_table
}
candidate_genes_for_cox_c3<-c(rownames(ph_hypo_table_c_3)[ph_hypo_table_c_3[,3]>0.05])
risk_score_table_multi_cox_c3<-riskscore(survival_cancer,
                                       candidate_genes_for_cox_c3,
                                       multi_variate_cox_c_3)

# make ROC curve
multi_ROC<-function(time_vector,risk_score_table){
  library(survivalROC)
  single_ROC<-function(single_time){
    for_ROC<-survivalROC(Stime = risk_score_table$DFS_time,
                         status = risk_score_table$DFS_status,
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
for_multi_ROC<-multi_ROC(time_vector = c(12*seq(3,5,0.2)),
                         risk_score_table = risk_score_table_multi_cox_c3)
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
high_low<-(risk_score_table_multi_cox_c3$total_risk_score>cut.off)
high_low[high_low==T]<-'high'
high_low[high_low==F]<-'low'
risk_score_table_multi_cox_c3<-cbind(risk_score_table_multi_cox_c3,high_low)

# prepare a dataset that AUC_max_time = 60
risk_score_table_multi_cox_c3_60<-risk_score_table_multi_cox_c3

# KM_plot generation
library(survminer)
# first edit the status of patients with OS > AUC max time. (censoring status=0 (alive), OS=365*5 days)
risk_score_table_multi_cox_c3$DFS_status[which(risk_score_table_multi_cox_c3$DFS_time>AUC_max_time)]<-0
risk_score_table_multi_cox_c3$DFS_time[which(risk_score_table_multi_cox_c3$DFS_time>AUC_max_time)]<-AUC_max_time
fit_km<-survfit(Surv(DFS_time,DFS_status)~high_low,data = risk_score_table_multi_cox_c3)
ggsurvplot(fit_km,conf.int = T,pval = T,legend.title='Total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),surv.median.line = 'hv',
           xlab="Time (month)",ylab="Cumulative DFS (percentage)",
           break.x.by=12)

ggsurvplot(fit_km, 
           conf.int = TRUE, # 增加置信区间
           pval = T,
           legend.title='Total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),
           fun = "cumhaz",
           xlab="Time (month)", # 绘制累计风险曲线
           break.x.by=12) 

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
                       AUC_max_time, ' months', sep = ''))
pROC
# AUCmax=0.71


# AUC_max_time=60
# KM_plot generation
library(survminer)
# first edit the status of patients with OS > AUC max time. (censoring status=0 (alive), OS=365*5 days)
risk_score_table_multi_cox_c3_60$DFS_status[which(risk_score_table_multi_cox_c3_60$DFS_time>60)]<-0
risk_score_table_multi_cox_c3_60$DFS_time[which(risk_score_table_multi_cox_c3_60$DFS_time>60)]<-60
fit_km_60<-survfit(Surv(DFS_time,DFS_status)~high_low,data = risk_score_table_multi_cox_c3_60)
ggsurvplot(fit_km_60,conf.int = T,pval = T,legend.title='Total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),surv.median.line = 'hv',
           xlab="Time (month)",ylab="Cumulative DFS (percentage)")

ggsurvplot(fit_km_60, 
           conf.int = TRUE, # 增加置信区间
           pval = T,
           legend.title='Total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),
           fun = "cumhaz",
           xlab="Time (month)") # 绘制累计风险曲线

# visualization of the ROC curves of multiple time point
pROC<-ggplot(for_multi_ROC,
             aes(x=False_positive,y=True_positive,label=Cut_values,color=Time_point))+
  geom_roc(labels = F,stat = 'identity',n.cuts = 0)+
  geom_abline(slope = 1,intercept = 0,color='red',linetype=2)+
  theme_bw()+
  theme(panel.background = linewidth(colour = 'black',size=1,fill='white'),
        panel.grid = element_blank())+
  annotate("text",x=0.75,y=0.15,
           label=paste("AUC max = ",round(AUC_max,2),'\n', 'AUC max time =', 
                       AUC_max_time, ' months', sep = ''))
pROC
# AUCmax=0.71


####alternative lasso optimization####
# alternative starting fromt the cross validation
# default type.measure="deviance", here try C-index.Both are good for cox model.
fitcv_test <- cv.glmnet(X,Y,family="cox", alpha=1,nfolds=10,type.measure = 'C')
plot(fitcv_test)
coef(fitcv_test, s="lambda.min")
best_lambda_test<-fitcv_test$lambda.min
best_lambda_test
# result in  14 or 10 variables, choose the result of deviance instead

# test the best model
best_model_test<-glmnet(X,Y,family='cox',alpha=1,lambda = best_lambda_test)
coef(best_model_test)
print(best_model_test)

# demonstrate the cox model
coxm_test<-cph(Surv(DFS_time,DFS_status==1)~FN1+LGALS3+SFTPB+HSPA1A+SERPINA1+PLCG2+HSPD1+MAPK13+HMGB2+STMN1+PASK+PPP1R7,
               x=T,y=T,data=survival_cancer,surv=T)
cox.zph(coxm_test)
#         chisq   df  p
#FN1      2.92081  1 0.087
#LGALS3   2.11814  1 0.146
#SFTPB    1.91871  1 0.166
#HSPA1A   0.57109  1 0.450
#SERPINA1 2.16657  1 0.141
#PLCG2    0.12709  1 0.721
#HSPD1    0.10253  1 0.749
#MAPK13   0.41711  1 0.518
#HMGB2    0.92183  1 0.337
#STMN1    0.17469  1 0.676
#PASK     0.00967  1 0.922
#PPP1R7   0.18720  1 0.665
#GLOBAL   8.83786 12 0.717
coxm_test
plot(cox.zph(coxm_test))
ggcoxzph(cox.zph(coxm_test))
# ideally the cox plot should be a line when k=0, suggesting that y won't change as time changes
ggcoxdiagnostics(coxm_test, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# Specifying the argument type = “dfbeta”, plots the estimated changes in the regression coefficients upon deleting each observation in turn
ggcoxdiagnostics(coxm_test, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
# 'The pattern looks fairly symmetric around 0' is the best result.
# FN1+LGALS3+SFTPB+HSPA1A+SERPINA1+PLCG2+HSPD1+MAPK13+HMGB2+STMN1+PASK+PPP1R7 are identified

ggpairs(survival_cancer[,c('FN1','LGALS3','SFTPB','HSPA1A','SERPINA1','PLCG2','HSPD1','MAPK13','HMGB2','STMN1','PASK','PPP1R7')],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())

formula_for_multivariate_ctest<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                              paste(c('FN1','LGALS3','SFTPB','HSPA1A','SERPINA1','PLCG2','HSPD1','MAPK13','HMGB2','STMN1','PASK','PPP1R7'),
                                                    sep='',collapse = '+')))
multi_variate_cox_c_3test<-coxph(formula_for_multivariate_ctest,data = survival_cancer)
ph_hypo_multi_c_3test<-cox.zph(multi_variate_cox_c_3test)
ph_hypo_table_c_3test<-ph_hypo_multi_c_3test$table[-nrow(ph_hypo_multi_c_3test$table),]

# forest plot
# HR<1 for protection, HR>1 for risk
ggforest(model = multi_variate_cox_c_3test,
         data = survival_cancer,
         main = "Hazard ratios of candidate genes",
         fontsize = 1)

# evaluate the cox model
## concordance index
C_indextest<-multi_variate_cox_c_3test$concordance['concordance']
if(C_indextest>=0.9){
  print('High accuracy')
}else{
  if(C_indextest<0.9&C_index>=0.7){
    print('Medium accuracy')
  }else{
    print('Low accuracy')
  }
}
C_indextest
# "Medium accuracy" 0.83

candidate_genes_for_cox_c3test<-c(rownames(ph_hypo_table_c_3test)[ph_hypo_table_c_3test[,3]>0.05])
risk_score_table_multi_cox_c3test<-riskscore(survival_cancer,
                                         candidate_genes_for_cox_c3test,
                                         multi_variate_cox_c_3test)

# make ROC curve
multi_ROC<-function(time_vector,risk_score_table){
  library(survivalROC)
  single_ROC<-function(single_time){
    for_ROC<-survivalROC(Stime = risk_score_table$DFS_time,
                         status = risk_score_table$DFS_status,
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
for_multi_ROCtest<-multi_ROC(time_vector = c(12*seq(3,5,0.2)),
                         risk_score_table = risk_score_table_multi_cox_c3test)
# multiple lines are overlapped

# KM-curve
# use the time with max AUC to do KM analysis
AUC_maxtest<-max(for_multi_ROCtest$AUC)
# select the last time point with max AUC indicating longer survival
AUC_max_timetest<-for_multi_ROCtest$Time_point[which(for_multi_ROCtest$AUC == AUC_maxtest)]
# remove duplicate
AUC_max_timetest<-AUC_max_timetest[!duplicated(AUC_max_timetest)]
AUC_max_timetest<-AUC_max_timetest[length(AUC_max_timetest)]
for_multi_ROCtest$Time_point<-as.factor(for_multi_ROCtest$Time_point)
# find the optimal cutoff value within the ROC curve of the optimal time point
optimal_time_ROC_dftest<-for_multi_ROCtest[which(for_multi_ROCtest$Time_point==AUC_max_timetest),]
cut.offtest<-optimal_time_ROC_dftest$Cut_values[which.max(optimal_time_ROC_dftest$True_positive-optimal_time_ROC_dftest$False_positive)]
high_lowtest<-(risk_score_table_multi_cox_c3test$total_risk_score>cut.offtest)
high_lowtest[high_lowtest==T]<-'high'
high_lowtest[high_lowtest==F]<-'low'
risk_score_table_multi_cox_c3test<-cbind(risk_score_table_multi_cox_c3test,high_lowtest)

# KM_plot generation
library(survminer)
# first edit the status of patients with OS > AUC max time. (censoring status=0 (alive), OS=365*5 days)
risk_score_table_multi_cox_c3test$DFS_status[which(risk_score_table_multi_cox_c3test$DFS_time>AUC_max_timetest)]<-0
risk_score_table_multi_cox_c3test$DFS_time[which(risk_score_table_multi_cox_c3test$DFS_time>AUC_max_timetest)]<-AUC_max_timetest
fit_kmtest<-survfit(Surv(DFS_time,DFS_status)~high_lowtest,data = risk_score_table_multi_cox_c3test)
ggsurvplot(fit_kmtest,conf.int = T,pval = T,legend.title='Total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),surv.median.line = 'hv',
           xlab="Time (month)",ylab="Cumulative DFS (percentage)")

ggsurvplot(fit_kmtest, 
           conf.int = TRUE, # 增加置信区间
           pval = T,
           legend.title='Total risk score',
           legend.labs=c('high','low'),risk.table = T,
           palette = c('red','blue'),
           fun = "cumhaz",
           xlab="Time (month)") # 绘制累计风险曲线

# visualization of the ROC curves of multiple time point
pROCtest<-ggplot(for_multi_ROCtest,
             aes(x=False_positive,y=True_positive,label=Cut_values,color=Time_point))+
  geom_roc(labels = F,stat = 'identity',n.cuts = 0)+
  geom_abline(slope = 1,intercept = 0,color='red',linetype=2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill='white'),
        panel.grid = element_blank())+
  annotate("text",x=0.75,y=0.15,
           label=paste("AUC max = ",round(AUC_maxtest,2),'\n', 'AUC max time =', 
                       AUC_max_timetest, ' months', sep = ''))
pROCtest

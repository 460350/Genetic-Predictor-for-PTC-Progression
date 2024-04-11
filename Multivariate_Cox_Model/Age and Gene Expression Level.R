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
library(ggsignif)
library(ggpubr)

# To evaluate if MAPK13, STMN1, and HMGB2 expression are associated with patients' age
# Prepare the clinical metadata for age stratification
summary(clin$`Diagnosis Age`)
# Min.    1st Qu. Median  Mean    3rd Qu.  Max. 
# 15.00   35.00   46.00   47.28   58.00   89.00 
clin_meta_age<-clin[,c(
  'Sample ID',
  'Diagnosis Age',
  'Neoplasm Disease Stage American Joint Committee on Cancer Code',
  'Cancer Type Detailed',
  'Disease Free (Months)',
  'Disease Free Status',
  'American Joint Committee on Cancer Metastasis Stage Code',
  'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code',
  'American Joint Committee on Cancer Tumor Stage Code',
  'Person Neoplasm Cancer Status',
  'Primary Lymph Node Presentation Assessment'
)]
dim(clin_meta_age)
rownames(clin_meta_age)<-clin_meta_age$`Sample ID`
# rename colname
colnames(clin_meta_age)<-c('ID','Diagnosis_age','Stage','Subtype','DFS_time','DFS_status',
                       'M_stage','N_stage','T_stage','Tumor_status','Primary_LN_status')

# define missing value as NA
clin_meta_age[clin_meta_age==""]<-NA

table(str_sub(rownames(clin_meta_age),14,15))
# 01    06 
# 499   1 

# remove the sample with ID ends with 06
n<-str_sub(rownames(clin_meta_age),15,15)!=6
clin_meta_age<-clin_meta_age[n,]
ncol(clin_meta_age)
table(str_sub(rownames(clin_meta_age),14,15))
#  01 
# 499

# modified the stage information
library(stringr)
head(clin_meta_age$Stage)
table(clin_meta_age$Stage)
# STAGE I  STAGE II STAGE III  STAGE IV STAGE IVA STAGE IVC 
#  282        51       110         2        46         6 

a = str_extract_all(clin_meta_age$Stage,"I|V");head(a)
b = sapply(a,paste,collapse = "");head(b)

clin_meta_age$simp_stage<-b
table(clin_meta_age$simp_stage)
#  I   II III  IV  NA 
# 282  51 110  54   2 

# remove the samples with NA stage info
clin_meta_age<-clin_meta_age[clin_meta_age$simp_stage!='NA',]

# make the clinical metadata based on DFS
clin_meta_age<-clin_meta_age[clin_meta_age$DFS_status!='NA',]
# Match the clin_data and exp
s = intersect(rownames(clin_meta_age),colnames(exp));length(s)
exp_age<-exp[,s]
colnames(clin_meta_age)
clin_meta_age<-clin_meta_age[s,]
dim(clin_meta_age)

# change the DFS_status
clin_meta_age$DFS_status<-ifelse(clin_meta_age$DFS_status=='0:DiseaseFree',
                                 0,
                                 1)
table(clin_meta_age$DFS_status)
#  0   1 
# 327  25 

# remove the samples without survival time
k1<-clin_meta_age$DFS_time>0;table(k1)
# FALSE  TRUE 
#  1     351
clin_meta_age<-clin_meta_age[k1,]
exp_age<-exp_age[,k1]

# make the exp matrix to cpm
exp_age<-log2(edgeR::cpm(as.matrix(exp_age))+1)


# evaluate the age distribution in TCGA model
summary(clin_meta_age$Diagnosis_age)
ggplot(data=clin_meta_age[,c(1,2)],aes(x=Diagnosis_age))+geom_histogram()

# group the data based on age ≤30. Age ≤30 = 1, otherwise 0.
clin_meta_age$age30<-ifelse(clin_meta_age$Diagnosis_age<=30,1,0)
table(clin_meta_age$age30)
#  0   1 
# 288  63

# prepare the survival table with patient's age and age30
survival_cancer_age<-cbind(survival_cancer,clin_meta_age[,c(2,13)])
survival_cancer_age$age30group<-ifelse(clin_meta_age$Diagnosis_age<=30,'≤30','>30')
survival_cancer_3gene<-survival_cancer_age[,c('TCGA_IDs','Diagnosis_age','age30','age30group',
                                            'DFS_time','DFS_status','MAPK13','STMN1','HMGB2')]

# evaluate the MAPK13, STMN1, and HMGB2 expression in patients age ≤30 and >30
# MAPK13
p1<-ggplot(survival_cancer_3gene,aes(x=age30group,y=MAPK13,fill=age30group))+
  geom_boxplot(aes(group=age30group),show.legend = T)+
  labs(fill='Age',x="Age")+
  scale_fill_hue(labels=c('Age>30','Age≤30'))+
  scale_fill_manual(values=c('red','blue'))+
  theme_classic()
p1

t.testResult1<-t.test(survival_cancer_3gene$MAPK13[survival_cancer_3gene$age30==0],
                      survival_cancer_3gene$MAPK13[survival_cancer_3gene$age30==1],
                      var.equal=T)
t.testResult1$p.value
# 0.0006751049

# STMN1
p2<-ggplot(survival_cancer_3gene,aes(x=age30group,y=STMN1,fill=age30group))+
  geom_boxplot(show.legend = T)+
  labs(fill='Age',x="Age")+
  scale_fill_hue(labels=c('Age>30','Age≤30'))+
  scale_fill_manual(values=c('red','blue'))+
  theme_classic()
p2
t.testResult2<-t.test(survival_cancer_3gene$STMN1[survival_cancer_3gene$age30==0],
                      survival_cancer_3gene$STMN1[survival_cancer_3gene$age30==1],
                      var.equal=T)
t.testResult2$p.value
# 0.6245512
# HMGB2
p3<-ggplot(survival_cancer_3gene,aes(x=age30group,y=HMGB2,fill=age30group))+
  geom_boxplot(show.legend = T)+
  labs(fill='Age',x='Age')+
  scale_fill_hue(labels=c('Age>30','Age≤30'))+
  scale_fill_manual(values=c('red','blue'))+
  theme_classic()
p3
t.testResult3<-t.test(survival_cancer_3gene$HMGB2[survival_cancer_3gene$age30==0],
                      survival_cancer_3gene$HMGB2[survival_cancer_3gene$age30==1],
                      var.equal=T)
t.testResult3$p.value
# 0.4149742

# combine 3 figures together
ggarrange(p1,p2,p3,ncol = 3,nrow = 1)

# evaluate if age30 could predict DFS in TCGA-THCA
t.testResult4<-t.test(survival_cancer_3gene$DFS_status[survival_cancer_3gene$age30==0],
                      survival_cancer_3gene$DFS_status[survival_cancer_3gene$age30==1],
                      var.equal=T)
t.testResult4$p.value
# 0.7822826

# evaluate if age difference between DFS=1 and DFS=0
t.testResult6<-t.test(survival_cancer_3gene$Diagnosis_age[survival_cancer_3gene$DFS_status==0],
                      survival_cancer_3gene$Diagnosis_age[survival_cancer_3gene$DFS_status==1],
                      var.equal=T)
t.testResult6$p.value
# 0.6370566

#### age difference between 3 gene high and 3 gene low ####
# add the 3 gene info from risk_score_table_multi_cox_c3
survival_cancer_3gene<-cbind(survival_cancer_3gene,risk_score_table_multi_cox_c3[,8])
colnames(survival_cancer_3gene)[10]<-'GeneStatus'
# t.test between age
p5<-ggplot(survival_cancer_3gene,aes(x=GeneStatus,y=Diagnosis_age,fill=GeneStatus))+
  geom_boxplot(aes(group=GeneStatus),show.legend = T)+
  scale_fill_hue(labels=c('High','Low'))+
  scale_fill_manual(values=c('red','blue'))+
  theme_classic()
p5

t.testResult5<-t.test(survival_cancer_3gene$Diagnosis_age[survival_cancer_3gene$GeneStatus=='low'],
                      survival_cancer_3gene$Diagnosis_age[survival_cancer_3gene$GeneStatus=='high'],
                      var.equal=T)
t.testResult5$p.value
# 0.09728524

####using age for multivariate analysis####
formula_for_multivariate_age<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                              paste(gsub(c(C_list,'Diagnosis_age'),
                                                         pattern = '-',replacement='_'),sep='',collapse = '+')))
multi_variate_cox_age<-coxph(formula_for_multivariate_age,data = survival_cancer_age)
# check if variances are supported by PH hypothesis
ph_hypo_multi_age<-cox.zph(multi_variate_cox_age)
# The last row of the table records the test results on the GLOBAL model. Delete it.
ph_hypo_table_age<-ph_hypo_multi_age$table[-nrow(ph_hypo_multi_age$table),]
# remove variances not supported by ph hypothesis and perform the 2nd regression
formula_for_multivariate_age<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                              paste(rownames(ph_hypo_table_age)[ph_hypo_table_age[,3]>0.05],
                                                    sep='',collapse = '+')))
multi_variate_cox_age_2<-coxph(formula_for_multivariate_age,data = survival_cancer_age)
ph_hypo_multi_age_2<-cox.zph(multi_variate_cox_age_2)
ph_hypo_table_age_2<-ph_hypo_multi_age_2$table[-nrow(ph_hypo_multi_age_2$table),]

# check the co-linearity between samples
correlation_age<-cor(survival_cancer_age[,rownames(ph_hypo_table_age_2)[ph_hypo_table_age_2[,3]>0.05]],
                   method = 'pearson')
library(GGally)
ggpairs(survival_cancer_age[,rownames(ph_hypo_table_age_2)[ph_hypo_table_age_2[,3]>0.05]],
        axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',size=1,fill = 'white'),
        panel.grid = element_blank())
library(rms)
vif_age<-rms::vif(multi_variate_cox_age_2)
sqrt(vif_age)<3

# some variable the colinearirity is significant, reduce it with lasso
# define response variable
Y<-as.matrix(survival::Surv(survival_cancer_age$DFS_time,survival_cancer_age$DFS_status))

# define a matrix of predictor variables
X<-as.matrix(survival_cancer_age[,gsub(rownames(ph_hypo_table_age_2),
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

# test the best model
best_model<-glmnet(X,Y,family='cox',alpha=1,lambda = best_lambda)
coef(best_model)
print(best_model)

# demonstrate the cox model
coxm_age<-cph(Surv(DFS_time,DFS_status==1)~MAPK13+HMGB2+STMN1+NMU,x=T,y=T,data=survival_cancer_age,surv=T)
cox.zph(coxm_age)
# chisq    df       p
# MAPK13 0.902  1 0.34
# HMGB2  2.071  1 0.15
# STMN1  0.625  1 0.43
# NMU    4.211  1 0.04
# GLOBAL 6.422  4 0.17
coxm_age

formula_for_multivariate_age<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                              paste(c('MAPK13','HMGB2','STMN1',"NMU"),
                                                    sep='',collapse = '+')))
multi_variate_cox_age_3<-coxph(formula_for_multivariate_age,data = survival_cancer_age)
ph_hypo_multi_age_3<-cox.zph(multi_variate_cox_age_3)
ph_hypo_table_age_3<-ph_hypo_multi_age_3$table[-nrow(ph_hypo_multi_age_3$table),]
plot(cox.zph(multi_variate_cox_age_3))
ggcoxzph(cox.zph(multi_variate_cox_age_3))

# forest plot
# HR<1 for protection, HR>1 for risk
ggforest(model = multi_variate_cox_age_3,
         data = survival_cancer_age,
         main = "Hazard ratios of candidate genes",
         fontsize = 1)

# add diagnosis age with MAPK13, HMGB2, and STMN1 and do the Cox regression analysis to draw a forest plot
formula_for_multivariate_age<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                                paste(c('MAPK13','HMGB2','STMN1',"Diagnosis_age"),
                                                      sep='',collapse = '+')))
multi_variate_cox_age_4<-coxph(formula_for_multivariate_age,data = survival_cancer_age)
ph_hypo_multi_age_4<-cox.zph(multi_variate_cox_age_4)
ph_hypo_table_age_4<-ph_hypo_multi_age_4$table[-nrow(ph_hypo_multi_age_4$table),]
plot(cox.zph(multi_variate_cox_age_4))
ggcoxzph(cox.zph(multi_variate_cox_age_4))

# forest plot
# HR<1 for protection, HR>1 for risk
ggforest(model = multi_variate_cox_age_4,
         data = survival_cancer_age,
         main = "Hazard ratios of candidate factors",
         fontsize = 1)
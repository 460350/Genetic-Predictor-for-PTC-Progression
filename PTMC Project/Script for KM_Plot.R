library(vroom)
library(tidyverse)
library(survival)
library(survminer)

# import the data
PTMCData<-vroom('PTMC_FollowUp_DataSheet_Reduction.csv')
PTMCData<-as.data.frame(PTMCData)
str(PTMCData)
# further modified the data
colnames(PTMCData)[colnames(PTMCData) == "UnifocalinOneLobe"] <- "MultifocalinOneLobe"
PTMCData$MultifocalinBothLobes[is.na(PTMCData$MultifocalinBothLobes)]<-0
PTMCData$MultifocalinOneLobe[is.na(PTMCData$MultifocalinOneLobe)]<-0
PTMCData$Multifocality<-PTMCData$MultifocalinBothLobes+PTMCData$MultifocalinOneLobe
PTMCData$Surgery[is.na(PTMCData$Surgery)]<-0
PTMCData$`Progression(withoutTwoSmaller)`[is.na(PTMCData$`Progression(withoutTwoSmaller)`)]<-0
PTMCData$PregnancywithTumor[is.na(PTMCData$PregnancywithTumor)]<-0
PTMCData$`Progression(deltaD≥3mmorLN+)`[is.na(PTMCData$`Progression(deltaD≥3mmorLN+)`)]<-0
PTMCData$HashimotoThyroiditis[is.na(PTMCData$HashimotoThyroiditis)]<-0
PTMCData$ClosetoMembrane[is.na(PTMCData$ClosetoMembrane)]<-0

####KMploting####
# KMplot for surgery
ptmcSurgery<-PTMCData[,c(1,3,6:8,11:13,22)]
colnames(ptmcSurgery)[colnames(ptmcSurgery) == "TimetoSurgeryorLastFollowUp"] <- "Time"
colnames(ptmcSurgery)[colnames(ptmcSurgery) == "Age of Diagnosis"] <- "AgeofDiagnosis"
colnames(ptmcSurgery)[colnames(ptmcSurgery) == "Tumor Size at Diagnosis"] <- "TumorSizeatDiagnosis"
ptmcSurgery$Age55 <- ifelse(ptmcSurgery$`AgeofDiagnosis` > 55, 1, 0)
ptmcSurgery$oneRiskFactor<-ifelse(rowSums(ptmcSurgery[,c(5:7,9)])>0,1,0)
ptmcSurgery$TwoRiskFactor<-ifelse(rowSums(ptmcSurgery[,c(5:7,9)])>1,1,0)
str(ptmcSurgery)
# attachment
attach(ptmcSurgery)
Surv_Obj<-Surv(Time,Surgery)
fit<-survfit(Surv_Obj ~ Age55,
             data = ptmcSurgery)
print(fit)
summary(fit)
ggsurvplot(fit,data = ptmcSurgery,
           ylab="nonSurgery probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.x.by=25,
           legend.labs=c("Age<=55","Age>55"),
           palette = "lancet",
           pval = T,
           risk.table = T)

ggsurvplot(fit,data = ptmcSurgery,
           ylab="Surgery probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.x.by=25,
           legend.labs=c("Age<=55","Age>55"),
           fun="cumhaz",
           pval = T,
           risk.table = T,
           palette = "lancet")

fit1<-survfit(Surv_Obj ~ 1,
              data = ptmcSurgery)
print(fit1)
summary(fit1)
ggsurvplot(fit1,data = ptmcSurgery,
           ylab="nonSurgery probability",
           xlab="Time (month)",
           risk.table = T,
           palette = "lancet")
ggsurvplot(fit1,data = ptmcSurgery,
           ylab="Surgery probability",
           xlab="Time (month)",
           fun="cumhaz",
           conf.int = F,
           break.x.by=25,
           ylim=c(0,1),
           xlim=c(0,100),
           risk.table = T,
           palette = "lancet")

fit5<-survfit(Surv_Obj ~ oneRiskFactor,
              data = ptmcSurgery)
ggsurvplot(fit5,data = ptmcSurgery,
           ylab="nonSurgery probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("0 risk factors","≥1 risk factor"),
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit5,data = ptmcSurgery,
           ylab="Surgery probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("0 risk factors","≥1 risk factor"),
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)

fit7<-survfit(Surv_Obj ~ TwoRiskFactor,
              data = ptmcSurgery)
ggsurvplot(fit7,data = ptmcSurgery,
           ylab="nonSurgery probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("0/1 risk factors","≥2 risk factor"),
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit7,data = ptmcSurgery,
           ylab="Surgery probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("0/1 risk factors","≥2 risk factor"),
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)
# KMplot for progression
ptmcProgression<-PTMCData[,c(1:4,6:8,11:13,20,22,24)]
colnames(ptmcProgression)[colnames(ptmcProgression) == "TimetoSurgeryorLastFollowUp"] <- "Time"
colnames(ptmcProgression)[colnames(ptmcProgression) == "Progression(withoutTwoSmaller)"] <- "Progression"
colnames(ptmcProgression)[colnames(ptmcProgression) == "Age of Diagnosis"] <- "AgeofDiagnosis"
colnames(ptmcProgression)[colnames(ptmcProgression) == "Tumor Size at Diagnosis"] <- "TumorSizeatDiagnosis"
str(ptmcProgression)
ptmcProgression$Age55 <- ifelse(ptmcProgression$AgeofDiagnosis > 55, 1, 0)
ptmcProgression$oneRiskFactor<-ifelse(rowSums(ptmcProgression[,c(7:9,13,14)])>0,1,0)
ptmcProgression$TwoRiskFactor<-ifelse(rowSums(ptmcProgression[,c(7:9,13,14)])>1,1,0)
ptmcProgression$Age30 <- ifelse(ptmcProgression$AgeofDiagnosis > 30, 0, 1)
ptmcProgression$Size6 <- ifelse(ptmcProgression$TumorSizeAtDiagnosis > 6, 1, 0)
ptmcProgression$Gender <- ifelse(ptmcProgression$Gender == "F",1,0)
# attachment
attach(ptmcProgression)
Surv_Obj2<-Surv(Time,Progression)
fit2<-survfit(Surv_Obj2 ~ 1,
             data = ptmcProgression)
print(fit2)
summary(fit2)
ggsurvplot(fit2,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           risk.table = T)
ggsurvplot(fit2,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.25,
           break.x.by=10,
           conf.int = F,
           fun="cumhaz",
           risk.table = T,
           palette = "lancet")

# divided by age=55
fit3<-survfit(Surv_Obj2 ~ Age55,
              data = ptmcProgression)
ggsurvplot(fit3,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("Age<=55","Age>55"),
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit3,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("Age<=55","Age>55"),
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)

# divided by age=30
fit8<-survfit(Surv_Obj2 ~ Age30,
              data = ptmcProgression)
ggsurvplot(fit8,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("Age>30","Age≤30"),
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit8,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("Age>30","Age≤30"),
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)

# divided by TumorSizeatDiagnosis=6 mm
fit11<-survfit(Surv_Obj2 ~ Size6,
              data = ptmcProgression)
ggsurvplot(fit11,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("Size≤6","Size>6"),
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit11,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("Size≤6","Size>6"),
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)

####Percentage of number of risk factors in patients####
PTMCData$rfPercentage<-PTMCData$HashimotoThyroiditis+PTMCData$ClosetoMembrane+
  PTMCData$PregnancywithTumor+PTMCData$Multifocality
summary(PTMCData$rfPercentage)
table(PTMCData$rfPercentage)

ptmcProgression$rfPercentage<-ptmcProgression$HashimotoThyroiditis+ptmcProgression$ClosetoMembrane+
  ptmcProgression$PregnancywithTumor+ptmcProgression$Multifocality
summary(ptmcProgression$rfPercentage)
table(ptmcProgression$rfPercentage)
############

# divided by patients with at least one risk factor and without any risk factor
fit4<-survfit(Surv_Obj2 ~ oneRiskFactor,
              data = ptmcProgression)
ggsurvplot(fit4,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit4,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)

# divided by patients with at least two risk factor and without any risk factor
fit6<-survfit(Surv_Obj2 ~ TwoRiskFactor,
              data = ptmcProgression)
ggsurvplot(fit6,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("0/1 risk factors","≥2 risk factor"),
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit6,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           legend.labs=c("0/1 risk factors","≥2 risk factor"),
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)
####Identification of factors association with Progression/Surgery####
# Progression
attach(ptmcProgression)
## univariate analysis
# try one at a time
sing.cox<-coxph(Surv_Obj2~TumorSizeAtDiagnosis,data=ptmcProgression)
sing.cox
summary(sing.cox)
# multiple at a time
covariates<-c("Age30","Size6","Gender","HashimotoThyroiditis","ClosetoMembrane","Multifocality")
univ_formulas<-sapply(covariates,
                      function(x) as.formula(paste('Surv_Obj2~',x)))
univ_models<-lapply(univ_formulas, function(x){coxph(x,data = ptmcProgression)})
# summarize the result into a dataframe
univ_results<-lapply(univ_models,
                     function(x){
                       x<-summary(x)
                       p.value<-signif(x$waldtest["pvalue"],digits=2)
                       wald.test<-signif(x$waldtest["test"],digits=2)
                       beta<-signif(x$coefficients[1],digits=2);#coefficient beta
                       HR<-signif(x$coefficients[2],digits=2);#exp(beta)
                       HR.confint.lower<-signif(x$conf.int[,"lower .95"],2)
                       HR.confint.upper<-signif(x$conf.int[,"upper .95"],2)
                       HR<-paste0(HR,"(",
                                  HR.confint.lower,"-",HR.confint.upper,")")
                       res<-c(beta,HR,wald.test,p.value)
                       names(res)<-c("beta","HR(95% CI for HR)","wald.test","p.value")
                       return(res)#return(exp(cbind(coef(x),confint(x))))
                     })
res<-t(as.data.frame(univ_results,check.names=F))
as.data.frame(res)
## multivariate analysis
pro.cox<-coxph(Surv_Obj2~Age30+
                 Size6+
                 Gender+
                 HashimotoThyroiditis+
                 ClosetoMembrane+
                 Multifocality,
               data = ptmcProgression)
pro.cox
summary(pro.cox)
ggsurvplot(survfit(pro.cox),data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           risk.table = T,
           palette = "lancet",
           pval = T)

# Surgery
attach(ptmcSurgery)
## univariate analysis
# try one at a time
sing.cox2<-coxph(Surv_Obj~TumorSizeAtDiagnosis,data=ptmcSurgery)
sing.cox2
summary(sing.cox2)
# multiple at a time
covariates<-c("AgeofDiagnosis","TumorSizeAtDiagnosis","HashimotoThyroiditis","ClosetoMembrane","Multifocality")
univ_formulas<-sapply(covariates,
                      function(x) as.formula(paste('Surv_Obj~',x)))
univ_models<-lapply(univ_formulas, function(x){coxph(x,data = ptmcSurgery)})
# summarize the result into a dataframe
univ_results<-lapply(univ_models,
                     function(x){
                       x<-summary(x)
                       p.value<-signif(x$waldtest["pvalue"],digits=2)
                       wald.test<-signif(x$waldtest["test"],digits=2)
                       beta<-signif(x$coefficients[1],digits=2);#coefficient beta
                       HR<-signif(x$coefficients[2],digits=2);#exp(beta)
                       HR.confint.lower<-signif(x$conf.int[,"lower .95"],2)
                       HR.confint.upper<-signif(x$conf.int[,"upper .95"],2)
                       HR<-paste0(HR,"(",
                                  HR.confint.lower,"-",HR.confint.upper,")")
                       res<-c(beta,HR,wald.test,p.value)
                       names(res)<-c("beta","HR(95% CI for HR)","wald.test","p.value")
                       return(res)#return(exp(cbind(coef(x),confint(x))))
                     })
res2<-t(as.data.frame(univ_results,check.names=F))
as.data.frame(res2)
## multivariate analysis
surg.cox<-coxph(Surv_Obj~AgeofDiagnosis+
                 TumorSizeAtDiagnosis+
                 HashimotoThyroiditis+
                 ClosetoMembrane+
                 Multifocality,
               data = ptmcSurgery)
surg.cox
summary(surg.cox)
ggsurvplot(survfit(surg.cox),data = ptmcSurgery,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           conf.int = F,
           risk.table = T,
           palette = "lancet",
           pval = T)


####Stratification of AgeofDiagnosis and TumorSizeatDiagnosis####
## TumorSizeatDiagnosis
# <5 is 0, [5,6) is 1, [6,7) is 2, [7,8) is 3, [8,9) is 4, >=9 is 5
ptmcProgression$SizeRank<-cut(ptmcProgression$TumorSizeAtDiagnosis, 
                             breaks=c(-Inf, 5, 6, 7, 8, 9, Inf), 
                             labels=c(0, 1, 2, 3, 4, 5))
fit9<-survfit(Surv_Obj2 ~ SizeRank,
              data = ptmcProgression)
ggsurvplot(fit9,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           legend.labs=c("Size≤5","5<Size≤6","6<Size≤7","7<Size≤8","8<Size≤9","Size>9"),
           conf.int = F,
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit9,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           legend.labs=c("Size≤5","5<Size≤6","6<Size≤7","7<Size≤8","8<Size≤9","Size>9"),
           conf.int = F,
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)
# alternatives to reduce the number of groups
# <5 is 0, [5,6) is 1, [6,7) is 2, >=7 is 3. So every group has more than 20 patients
ptmcProgression$SizeRank2<-cut(ptmcProgression$TumorSizeAtDiagnosis, 
                              breaks=c(-Inf, 5, 6, 7, Inf), 
                              labels=c(0, 1, 2, 3))
fit10<-survfit(Surv_Obj2 ~ SizeRank2,
              data = ptmcProgression)
ggsurvplot(fit10,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           legend.labs=c("Size≤5","5<Size≤6","6<Size≤7","Size>7"),
           conf.int = F,
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit10,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           legend.labs=c("Size≤5","5<Size≤6","6<Size≤7","Size>7"),
           conf.int = F,
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)
## AgeofDiagnosis
# =<30 is 0, (30,40] is 1, (40,50] is 2, >50 is 3
ptmcProgression$AgeRank<-cut(ptmcProgression$AgeofDiagnosis, 
                              breaks=c(-Inf, 30, 40, 50, Inf), 
                              labels=c(0, 1, 2, 3))
fit10<-survfit(Surv_Obj2 ~ AgeRank,
              data = ptmcProgression)
ggsurvplot(fit10,data = ptmcProgression,
           ylab="nonProgression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           legend.labs=c("Age≤30","30<Age≤40","40<Age≤50","Age>50"),
           conf.int = F,
           risk.table = T,
           palette = "lancet",
           pval = T)
ggsurvplot(fit10,data = ptmcProgression,
           ylab="Progression probability",
           xlab="Time (month)",
           xlim=c(0,100),
           ylim=c(0,1),
           break.y.by=0.1,
           break.x.by=20,
           legend.labs=c("Age≤30","30<Age≤40","40<Age≤50","Age>50"),
           conf.int = F,
           fun="cumhaz",
           risk.table = T,
           palette = "lancet",
           pval = T)
## plot for distribution of Diagnosis Age and Size at Diagnosis

####Female only Data####
# import the data
F_PTMCData<-vroom('F_Only_PTMC_FollowUp_DataSheet_Reduction.csv')
F_PTMCData<-as.data.frame(F_PTMCData)
str(F_PTMCData)
# further modified the data
colnames(F_PTMCData)[colnames(F_PTMCData) == "UnifocalinOneLobe"] <- "MultifocalinOneLobe"
F_PTMCData$MultifocalinBothLobes[is.na(F_PTMCData$MultifocalinBothLobes)]<-0
F_PTMCData$MultifocalinOneLobe[is.na(F_PTMCData$MultifocalinOneLobe)]<-0
F_PTMCData$Multifocality<-F_PTMCData$MultifocalinBothLobes+F_PTMCData$MultifocalinOneLobe
F_PTMCData$Surgery[is.na(F_PTMCData$Surgery)]<-0
F_PTMCData$`Progression(withoutTwoSmaller)`[is.na(F_PTMCData$`Progression(withoutTwoSmaller)`)]<-0
F_PTMCData$PregnancywithTumor[is.na(F_PTMCData$PregnancywithTumor)]<-0
F_PTMCData$`Progression(deltaD≥3mmorLN+)`[is.na(F_PTMCData$`Progression(deltaD≥3mmorLN+)`)]<-0
F_PTMCData$HashimotoThyroiditis[is.na(F_PTMCData$HashimotoThyroiditis)]<-0
F_PTMCData$ClosetoMembrane[is.na(F_PTMCData$ClosetoMembrane)]<-0

# focus on progression
F_ptmcProgression<-F_PTMCData[,c(1,2,3,4,6:8,11:13,20,22,24)]
colnames(F_ptmcProgression)[colnames(F_ptmcProgression) == "TimetoSurgeryorLastFollowUp"] <- "Time"
colnames(F_ptmcProgression)[colnames(F_ptmcProgression) == "Progression(withoutTwoSmaller)"] <- "Progression"
colnames(F_ptmcProgression)[colnames(F_ptmcProgression) == "Age of Diagnosis"] <- "AgeofDiagnosis"
colnames(F_ptmcProgression)[colnames(F_ptmcProgression) == "Tumor Size at Diagnosis"] <- "TumorSizeatDiagnosis"
str(F_ptmcProgression)
F_ptmcProgression$Age30 <- ifelse(F_ptmcProgression$AgeofDiagnosis > 30, 0, 1)
F_ptmcProgression$Size6 <- ifelse(F_ptmcProgression$TumorSizeAtDiagnosis > 6, 1, 0)

attach(F_ptmcProgression)
Surv_Obj3<-Surv(Time,Progression)
# univariate cox regression model
F_covariates<-c("Age30","Size6","PregnancywithTumor","HashimotoThyroiditis","ClosetoMembrane","Multifocality")
F_univ_formulas<-sapply(F_covariates,
                      function(x) as.formula(paste('Surv_Obj3~',x)))
F_univ_models<-lapply(F_univ_formulas, function(x){coxph(x,data = F_ptmcProgression)})
F_univ_results<-lapply(F_univ_models,
                     function(x){
                       x<-summary(x)
                       p.value<-signif(x$waldtest["pvalue"],digits=2)
                       wald.test<-signif(x$waldtest["test"],digits=2)
                       beta<-signif(x$coefficients[1],digits=2);#coefficient beta
                       HR<-signif(x$coefficients[2],digits=2);#exp(beta)
                       HR.confint.lower<-signif(x$conf.int[,"lower .95"],2)
                       HR.confint.upper<-signif(x$conf.int[,"upper .95"],2)
                       HR<-paste0(HR,"(",
                                  HR.confint.lower,"-",HR.confint.upper,")")
                       res<-c(beta,HR,wald.test,p.value)
                       names(res)<-c("beta","HR(95% CI for HR)","wald.test","p.value")
                       return(res)#return(exp(cbind(coef(x),confint(x))))
                     })
F_res<-t(as.data.frame(F_univ_results,check.names=F))
as.data.frame(F_res)

library(vroom)
library(tidyverse)
library(survival)
library(survminer)

# import the data
PTMCData<-vroom('PTMC_FollowUp_DataSheet.csv')
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
PTMCData$HashimatoThyroiditis[is.na(PTMCData$HashimatoThyroiditis)]<-0
PTMCData$ClosetoMembrane[is.na(PTMCData$ClosetoMembrane)]<-0

####KMploting####
# KMplot for surgery
ptmcSurgery<-PTMCData[,c(1,3,6:8,11:13,22)]
colnames(ptmcSurgery)[colnames(ptmcSurgery) == "TimetoSurgeryorLastFollowUp"] <- "Time"
colnames(ptmcSurgery)[colnames(ptmcSurgery) == "Age of Diagnosis"] <- "AgeofDiagnosis"
colnames(ptmcSurgery)[colnames(ptmcSurgery) == "Tumor Size at Diagnosis"] <- "TumorSizeatDiagnosis"
ptmcSurgery$Age55 <- ifelse(ptmcSurgery$`Age of Diagnosis` > 55, 1, 0)
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
ptmcProgression<-PTMCData[,c(1,3,6:8,11:13,20,22)]
colnames(ptmcProgression)[colnames(ptmcProgression) == "TimetoSurgeryorLastFollowUp"] <- "Time"
colnames(ptmcProgression)[colnames(ptmcProgression) == "Progression(withoutTwoSmaller)"] <- "Progression"
colnames(ptmcProgression)[colnames(ptmcProgression) == "Age of Diagnosis"] <- "AgeofDiagnosis"
colnames(ptmcProgression)[colnames(ptmcProgression) == "Tumor Size at Diagnosis"] <- "TumorSizeatDiagnosis"
str(ptmcProgression)
ptmcProgression$Age55 <- ifelse(ptmcProgression$AgeofDiagnosis > 55, 1, 0)
ptmcProgression$oneRiskFactor<-ifelse(rowSums(ptmcProgression[,c(5:7,10)])>0,1,0)
ptmcProgression$TwoRiskFactor<-ifelse(rowSums(ptmcProgression[,c(5:7,10)])>1,1,0)
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
           break.x.by=25,
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
           legend.labs=c("0 risk factors","≥1 risk factor"),
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
           legend.labs=c("0 risk factors","≥1 risk factor"),
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
sing.cox<-coxph(Surv_Obj2~TumorSizeatDiagnosis,data=ptmcProgression)
sing.cox
summary(sing.cox)
# multiple at a time
covariates<-c("AgeofDiagnosis","TumorSizeatDiagnosis","HashimatoThyroiditis","ClosetoMembrane","PregnancywithTumor","Multifocality")
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
pro.cox<-coxph(Surv_Obj2~AgeofDiagnosis+
                 TumorSizeatDiagnosis+
                 HashimatoThyroiditis+
                 ClosetoMembrane+
                 PregnancywithTumor+
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
sing.cox2<-coxph(Surv_Obj~TumorSizeatDiagnosis,data=ptmcSurgery)
sing.cox2
summary(sing.cox2)
# multiple at a time
covariates<-c("AgeofDiagnosis","TumorSizeatDiagnosis","HashimatoThyroiditis","ClosetoMembrane","PregnancywithTumor","Multifocality")
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
                 TumorSizeatDiagnosis+
                 HashimatoThyroiditis+
                 ClosetoMembrane+
                 PregnancywithTumor+
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

####Percentage of number of risk factors in patients####
PTMCData$rfPercentage<-PTMCData$HashimatoThyroiditis+PTMCData$ClosetoMembrane+
  PTMCData$PregnancywithTumor+PTMCData$Multifocality
summary(PTMCData$rfPercentage)
table(PTMCData$rfPercentage)

ptmcProgression$rfPercentage<-ptmcProgression$HashimatoThyroiditis+ptmcProgression$ClosetoMembrane+
  ptmcProgression$PregnancywithTumor+ptmcProgression$Multifocality
summary(ptmcProgression$rfPercentage)
table(ptmcProgression$rfPercentage)

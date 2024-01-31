library(tidyverse)
library(vroom)

# load the DFS samples
load('DFS_THCA_PanCancer.RData')

# SURVIVAL ANALYSIS OF SINGLE GENE
# analysis the impact of single gene expression on the survival
genelist<-c('FN1','KRT19','SERPINA1','LGALS3','MET')
library(survival)
library(survminer)
splots = list()

for (i in 1:length(genelist)) {
  x = clin_meta_DFS
  g = genelist[i]
  x$gene = ifelse(exp_DFS[g,]>median(exp_DFS[g,]),'high','low')
  sfit1 = survfit(Surv(DFS_time,DFS_status)~gene, data = x)
  splots[[i]] = ggsurvplot(sfit1, pval = T, palette = c("red","blue"),
                           data = x, legend = c(0.8,0.8),
                           title = genelist[i],risk.table = T)
}
splots

####SURVIVAL ANALYSIS OF COX KM PLOT####
# try the DEGs from Scissor analysis
library(glmnet)
library(GGally)
library(rms)
library(survivalROC)
library(plotROC)
library(stringr)
# prepare the survival table
survival_cancer<-clin_meta_DFS%>%select(DFS_time,DFS_status,simp_stage)
# in package 'survival', 0 represents nothing happens, 1 represents Recurred/Progressed
colnames(survival_cancer)<-c("DFS_time","DFS_status","tumor_stage")
survival_cancer$DFS_time<-as.numeric(survival_cancer$DFS_time)
survival_cancer$DFS_status<-as.numeric(survival_cancer$DFS_status)

# remove the NA column
survival_cancer<-survival_cancer[!is.na(survival_cancer),]
rownames(survival_cancer)
# keep the rows with real data
survival_cancer<-survival_cancer[1:351,]
survival_cancer$TCGA_IDs<-rownames(survival_cancer)
rownames(survival_cancer)

# add the expression data to survival_cancer
survival_cancer<-cbind(survival_cancer,t(exp_DFS)[survival_cancer$TCGA_IDs,])
# change the gene names containing '-'
colnames(survival_cancer)<-gsub(colnames(survival_cancer),pattern='-',replacement='_')


Scissor.de.gene<-as.data.frame(vroom('Scissor_stage_v2_ptc_de_genes.csv'))
Scissor.de.gene<-Scissor.de.gene[,-1]
rownames(Scissor.de.gene)<-Scissor.de.gene$SYMBOL

# remove the genes that are not included in the data.frame survival cancer
Scissor.de.gene<-Scissor.de.gene%>%
  dplyr::filter(SYMBOL %in% colnames(survival_cancer[,-c(1:4)]))

# univariate cox regression
# filter potential useful sig genes using univariate cox regression
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
uni_cox_df<-uni_cox_in_bulk(gene_list = Scissor.de.gene$SYMBOL,survival_info_df = survival_cancer)
uni_cox_sig_genes<-uni_cox_df$uni_cox_sig_genes
uni_cox_sig_genes
# 17 genes identified

# use Lasso regression
x<-as.matrix(survival_cancer[,gsub(Scissor.de.gene$SYMBOL,
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
# 0 genes identified

# use the candidates identifed above to construct multivariate model
# perform the multi-variates cox regression using qualified genes
formula_for_multivariate<-as.formula(paste0('Surv(DFS_time,DFS_status)~',
                                            paste(uni_cox_sig_genes,sep='',collapse = '+')))
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
vif<5

# all true

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
# "Medium accuracy" 0.8113321

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
candidate_genes_for_cox2<-c(rownames(ph_hypo_table_2)[ph_hypo_table_2[,3]>0.05])
risk_score_table_multi_cox2<-riskscore(survival_cancer,
                                       candidate_genes_for_cox2,
                                       multi_variate_cox_2)

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
risk_score_table_multi_cox2$DFS_status[which(risk_score_table_multi_cox2$DFS_time>AUC_max_time)]<-0
risk_score_table_multi_cox2$DFS_time[which(risk_score_table_multi_cox2$DFS_time>AUC_max_time)]<-AUC_max_time
fit_km<-survfit(Surv(DFS_time,DFS_status)~high_low,data = risk_score_table_multi_cox2)
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
# AUCmax=0.29

####SCISSOR Analysis####
# the possible link between scRNAseq and bulkseq based on survival data
library(Scissor)
library(Seurat)
# modify the expr_matrix
bulk_dataset<-survival_cancer[,5:14559]
bulk_dataset<-t(as.matrix(bulk_dataset))
# modify the phenotype/survival data
phenotype<-clin_meta_DFS[,c(4,5)]
colnames(phenotype)<-c('time','status')
THCA_DFSstatus<-phenotype$status
names(THCA_DFSstatus)<-rownames(phenotype)
tags<-c('DiseaseFree','RecurredProgressed')

# overview of the DFS time
ggplot(data=phenotype,aes(phenotype$time))+
  geom_histogram(col='white',binwidth = 10)+
  labs(title = 'Histogram for DFS time in TCGA-THCA PanCancer Atlas 2018',x='DFS time')+
  theme_classic()

# necessary function: Scissor_change. Change the RNA as SCT
Scissor_change <- function(bulk_dataset, sc_dataset, phenotype, tag = NULL,
                           alpha = NULL, cutoff = 0.2, family = c("gaussian","binomial","cox"),
                           Save_file = "Scissor_inputs.RData", Load_file = NULL){
  library(Seurat)
  library(Matrix)
  library(preprocessCore)
  
  
  if (is.null(Load_file)){
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    if (class(sc_dataset) == "Seurat"){
      sc_exprs <- as.matrix(sc_dataset@assays$SCT@data)
      network  <- as.matrix(sc_dataset@graphs$SCT_snn)
    }else{
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", verbose = F)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = F)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, verbose = F)
      network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    diag(network) <- 0
    network[which(network != 0)] <- 1
    
    dataset0 <- cbind(bulk_dataset[common,], sc_exprs[common,])         # Dataset before quantile normalization.
    dataset1 <- normalize.quantiles(dataset0)                           # Dataset after  quantile normalization.
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    
    Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    if (quality_check[3] < 0.01){
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    if (family == "binomial"){
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)){
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }else{
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    if (family == "gaussian"){
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)){
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }else{
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", paste(tmp[1:(length(z)-1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    if (family == "cox"){
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2){
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      }else{
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)
  }else{
    load(Load_file)
  }
  
  if (is.null(alpha)){
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  }
  for (i in 1:length(alpha)){
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, nlambda = 100, nfolds = min(10,nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial"){
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X)+1)])
    }else{
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))
    
    if (percentage < cutoff){
      break
    }
    cat("\n")
  }
  print("|**************************************************|")
  
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = family),
              Coefs = Coefs,
              Scissor_pos = Cell1,
              Scissor_neg = Cell2))
}

# apply the Scissor analysis
infos1 <- Scissor_change(bulk_dataset, fhs_thyrocyte, THCA_DFSstatus, tag = tags,alpha = 0.3, 
                         family = "binomial", Save_file = 'Scissor_THCA_DFS.RData')
# Scissor identifies 726 Scissor+ cells associated with worse survival 
# and 1081 Scissor- cells associated with good survival

# visualizing the Scissor result in Seurat
Scissor_select <- rep(0, ncol(fhs_thyrocyte))
names(Scissor_select) <- colnames(fhs_thyrocyte)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
fhs_thyrocyte <- AddMetaData(fhs_thyrocyte, metadata = Scissor_select, col.name = "Scissor_DFS")
DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_DFS', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_DFS', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1),
        split.by = 'annotation', ncol = 3)

table(Idents(fhs_thyrocyte),fhs_thyrocyte$Scissor_DFS)

save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

Idents(fhs_thyrocyte)<-"Scissor_DFS"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="1",ident.2="2",
                                min.pct = 0.25,logfc.threshold = log(2),min.diff.pct = 0.25)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL

# export the gene list
write.csv(Scissor.de.markers,'Scissor_DFS_ptc_de_genes.csv')
Idents(fhs_thyrocyte)<-"annotation"

rm(fhs_thyrocyte)

## Reliability significance test (heavily time-consuming, 15 h for 50 times)
# setting
# This p-value is less than 0.05, indicating that these associations are reliable. 
numbers <- length(infos1$Scissor_pos) + length(infos1$Scissor_neg)
# load the infos1 result
load('Scissor_THCA_DFS.RData')
# in reality nfold should be ~100
result1 <- reliability.test(X, Y, network, alpha = 0.3, family = "binomial", cell_num = numbers, n = 10, nfold = 20)
# the result mostly is p=NA & NaN when alpha = 0.05 or 0.005 in infos1

# cell level evaluation
evaluate_summary <- evaluate.cell('Scissor_THCA_DFS.RData', infos1, FDR = 0.05, bootstrap_n = 100)
all(evaluate_summary$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary$`Correlation > 0`)) > 50)

# SAVE the result with Scissor result
save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

####visualization of Scissor+ and Scissor- cells
# identify the different-expressed genes between Scissor+ and Scissor- cells
# alpha = 0.05
Idents(fhs_thyrocyte)<-"Scissor_DFS"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="1",ident.2="2",
                                min.pct = 0.25,logfc.threshold = log(2),min.diff.pct = 0.25)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL

# export the gene list
write.csv(Scissor.de.markers,'Scissor_DFS_ptc_de_genes.csv')
Idents(fhs_thyrocyte)<-"annotation"


# visualization with volcano plot
Scissor.de.markers<-Scissor.de.markers%>%
  mutate(Expression=case_when(avg_log2FC>=1&p_val_adj<0.05~"Up-regulated",
                              avg_log2FC<=-1&p_val_adj<0.05~"Down-regulated",
                              TRUE~"Unchanged"))
p1<-ggplot(Scissor.de.markers,aes(x=avg_log2FC,y=-log(p_val_adj,10)))+geom_point(aes(color=Expression),size=1)+
  xlab(expression("log"[2]*"FC"))+ylab(expression("-log"[10]*"FDR"))+
  scale_color_manual(values=c("dodgerblue3","gray50","firebrick3"))+
  guides(color=guide_legend(override.aes = list(size=1.5)))
p1
# identify the top 5 genes
top<-5
top_Sci_genes<-bind_rows(
  Scissor.de.markers%>%
    filter(Expression=='Up-regulated')%>%
    arrange(p_val_adj,desc(abs(avg_log2FC)))%>%
    head(top),
  Scissor.de.markers%>%
    filter(Expression=='Down-regulated')%>%
    arrange(p_val_adj,desc(abs(avg_log2FC)))%>%
    head(top),
)
p2<-p1+geom_label(data=top_Sci_genes,mapping=aes(x=avg_log2FC,y=-log(p_val_adj,10),
                                                 label=SYMBOL),size=3)
p2

####GSEA####
library(Seurat)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)

### enrichment analysis
geneList<-Scissor.de.markers[,2]
names(geneList)<-as.character(Scissor.de.markers[,6])
geneList<-sort(geneList,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is significant

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is significant
# select the KEGG pathway only
emC2_KEGG<-emC2
emC2_KEGG@result<-emC2_KEGG@result[grep("KEGG_", emC2@result$Description),]
enrichplot::dotplot(emC2_KEGG,showCategory = 5,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
emC5BP<-GSEA(geneList,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is significant

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
emC5CC<-GSEA(geneList,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is significant

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is significant

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'Scissor_pos & Scissor_neg GSEA Result scRNA.RData')

## GSEA for bulk RNAseq

#### Prepare C2KEGG for the NES plot
emC2_KEGG_mod<-emC2_KEGG
rownames(emC2_KEGG@result)
emC2_KEGG_mod@result<-emC2_KEGG@result[-c(3,4,6,7,8,9,10,12,14,16,17,23,24,29),]
enrichplot::dotplot(emC2_KEGG_mod,showCategory = 5,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.26,0.7)
p5<-ggplot(emC2_KEGG_mod@result%>%arrange(desc(NES))%>%head(n=10),aes(reorder(Description,NES),NES))+
  geom_col(aes(fill=NES))+
  coord_flip()+
  labs(x="KEGG",y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA")
p5
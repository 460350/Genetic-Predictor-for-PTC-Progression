library(vroom)
library(DESeq2)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(stringr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)

# prepare the countData
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

# Match the clin_data and exp
s = intersect(rownames(clin_meta_DFS),colnames(exp));length(s)
DFS_exp<-exp[,s]

countData<-DFS_exp

# gene filtration
# it has been filtered

# remove the duplicate sample
# reorder the expression matrix based on the sample ID
countData<-countData[,sort(colnames(countData))]
# identified the duplicated sample based the sample ID from 1 to 15
l<-!duplicated(str_sub(colnames(countData),1,15));table(l)
# TRUE 
# 351 
countData<-countData[,l]
colnames(countData)=str_sub(colnames(countData),1,15)
ncol(countData)

# prepare the metadata
# Identify the IDs indicating normal sample
DFS_Ident<-clin_meta_DFS$DFS_status
DFS_Ident<-as.factor(DFS_Ident)
names(DFS_Ident)<-colnames(countData)
DFS_Ident<-as.data.frame(DFS_Ident)
DFS_Ident$DFS_Ident<-as.factor(DFS_Ident$DFS_Ident)

# PCA analysis
t_countData<-as.data.frame(t(countData))
t_countData$tissue=DFS_Ident$DFS_Ident
t_countData_pca<-PCA(t_countData[,-14556],graph=T)
fviz_screeplot(t_countData_pca,addlabels=T)
var<-get_pca_var(t_countData_pca)
fviz_pca_ind(t_countData_pca,
             mean.point=F,#去除分组的中心点
             label="none",
             habillage=as.factor(t_countData$tissue),#根据样本类型着色,should be a factor
             addEllipses = T)#添加边界线

# transform the countData into integer
DFS_countData<-round(countData,0)
DFS_countData<-as.integer(DFS_countData)

# construct a DESeqDataSet object
# design means the column including treatment (design factor). The design column has to be factor
dds<-DESeqDataSetFromMatrix(countData=DFS_countData,colData=DFS_Ident,
                            design=~DFS_Ident)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
# make a logical value keep for the filtering
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]

# set a factor level. Untreated as reference level
dds$DFS_Ident<-relevel(dds$DFS_Ident,ref="0")

# run DESeq
dds<-DESeq(dds)
res<-results(dds,alpha=0.01)

res

# explore results
summary(res)
# out of 14555 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 89, 0.61%
# LFC < 0 (down)     : 43, 0.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 11)

res_result<-as.data.frame(res)
res_result=res_result[order(res_result$pvalue),]
res_result$sample_ID<-rownames(res_result)

res_result<-res_result%>%
  mutate(Expression=case_when(log2FoldChange>=log(2)&padj<0.05~"Up-regulated",
                              log2FoldChange<=-log(2)&padj<0.05~"Down-regulated",
                              TRUE~"Unchanged"))

res_result<-res_result%>%
  mutate(Significance=case_when(abs(log2FoldChange)>=log(2)&padj<0.05&padj>0.01~"FDR 0.05",
                                abs(log2FoldChange)>=log(2)&padj<=0.01&padj>0.001~"FDR 0.01",
                                abs(log2FoldChange)>=log(2)&padj<=0.001~"FDR 0.001",
                                TRUE~"Unchanged"))

# counting the upregulated and downregulated genes
res_count<-res_result%>%
  dplyr::count(Expression,Significance)
res_count
#      Expression Significance     n
# 1 Down-regulated    FDR 0.001    19
# 2 Down-regulated     FDR 0.01    20
# 3 Down-regulated     FDR 0.05    40
# 4      Unchanged    Unchanged 14392
# 5   Up-regulated    FDR 0.001    42
# 6   Up-regulated     FDR 0.01    21
# 7   Up-regulated     FDR 0.05    21

# identify the top 10 genes
top<-10
top_genes<-bind_rows(
  res_result%>%
    filter(Expression=='Up-regulated')%>%
    arrange(padj,desc(abs(log2FoldChange)))%>%
    head(top),
  res_result%>%
    filter(Expression=='Down-regulated')%>%
    arrange(padj,desc(abs(log2FoldChange)))%>%
    head(top),
)

# select the top genes seq data
top_t_countData<-t_countData%>%
  select(all_of(rownames(top_genes)))

# volcano plot
plotMA(res)

p3<-ggplot(res_result,aes(x=log2FoldChange,y=-log(padj,10)))+geom_point(aes(color=Expression),size=2/5)+
  xlab(expression("log"[2]*"FC"))+ylab(expression("-log"[10]*"FDR"))+
  scale_color_manual(values=c("dodgerblue3","gray50","firebrick3"))+
  guides(color=guide_legend(override.aes = list(size=1.5)))
p3

p4<-p3+geom_label(data=top_genes,mapping=aes(x=log2FoldChange,y=-log(padj,10),
                                             label=sample_ID),size=2)
p4

# hierarchical clustering
data_scaled<-scale(x=top_t_countData,center=T,scale=T)

#计算距离
res_dist<-dist(x=data_scaled,method="euclidean")
# print distance matrix
output=as.matrix(res_dist)[1:40,1:40]

# Hierarchical Clustering
res_hc<-hclust(d=res_dist,method="complete")
fviz_dend(x=res_hc,cex=0.7,lwd=0.7,k=4,
          k_colors = c("jco"),rect=T,rect_border="jco",rect_fill = T)


# enrichment analysis
# prepare data
geneList<-res_result[,2]
names(geneList)<-as.character(res_result[,7])
geneList<-sort(geneList,decreasing = T)

# GSEA
H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 10,split=".sign")+facet_grid(~.sign)
# select the KEGG pathway only
emC2_KEGG<-emC2
emC2_KEGG@result<-emC2_KEGG@result[grep("KEGG_", emC2@result$Description),]
enrichplot::dotplot(emC2_KEGG,showCategory = 20,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
emC5BP<-GSEA(geneList,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
emC5CC<-GSEA(geneList,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'DFS_pos & DFS_neg GSEA Result bulk.RData')

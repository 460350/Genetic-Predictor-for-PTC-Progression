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
countData<-exp

# gene filtration
# only the genes with count > 10 in >50% genes will be kept
k = apply(countData,1, function(x){sum(x>10)>0.5*ncol(countData)});table(k)
# FALSE  TRUE 
# 11172 19765 
countData = countData[k,]
nrow(countData)

# remove the duplicate sample
# reorder the expression matrix based on the sample ID
countData<-countData[,sort(colnames(countData))]
# identified the duplicated sample based the sample ID from 1 to 15
l<-!duplicated(str_sub(colnames(countData),1,15));table(l)
# TRUE 
# 570 
countData<-countData[,l]
colnames(countData)=str_sub(colnames(countData),1,15)
ncol(countData)

# prepare the metadata
# identify the tumor and normal sample
table(str_sub(colnames(countData),14,15))
# 01    06  11 
# 503   8   59 
# remove the sample with ID ends with 06
m<-str_sub(colnames(countData),15,15)!=6
countData<-countData[,m]
ncol(countData)
table(str_sub(colnames(countData),14,15))
# 01   11 
# 503  59 
# Identify the IDs indicating normal sample
TN_Ident<-ifelse(as.numeric(str_sub(colnames(countData),14,15)) < 10,'tumor','normal')
TN_Ident<-as.factor(TN_Ident)
names(TN_Ident)<-colnames(countData)
TN_Ident<-as.data.frame(TN_Ident)

# PCA analysis
t_countData<-as.data.frame(t(countData))
t_countData$tissue=TN_Ident$TN_Ident
t_countData_pca<-PCA(t_countData[,-19766],graph=T)
fviz_screeplot(t_countData_pca,addlabels=T)
var<-get_pca_var(t_countData_pca)
fviz_pca_ind(t_countData_pca,
             mean.point=F,#去除分组的中心点
             label="none",
             habillage=as.factor(t_countData$tissue),#根据样本类型着色,should be a factor
             addEllipses = T)#添加边界线

# construct a DESeqDataSet object
# design means the column including treatment (design factor). The design column has to be factor
dds<-DESeqDataSetFromMatrix(countData=countData,colData=TN_Ident,
                            design=~TN_Ident)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
# make a logical value keep for the filtering
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]

# set a factor level. Untreated as reference level
dds$TN_Ident<-relevel(dds$TN_Ident,ref="normal")

# run DESeq
dds<-DESeq(dds)
res<-results(dds,alpha=0.01)

res 

# explore results
summary(res)
# out of 19765 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 6550, 33%
# LFC < 0 (down)     : 6023, 30%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%

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
#       Expression Significance     n
# 1 Down-regulated    FDR 0.001  1898
# 2 Down-regulated     FDR 0.01    45
# 3 Down-regulated     FDR 0.05     6
# 4      Unchanged    Unchanged 15192
# 5   Up-regulated    FDR 0.001  2594
# 6   Up-regulated     FDR 0.01    24
# 7   Up-regulated     FDR 0.05     6

# identify the top 25 genes
top<-50
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

p1<-ggplot(res_result,aes(x=log2FoldChange,y=-log(padj,10)))+geom_point(aes(color=Expression),size=2/5)+
  xlab(expression("log"[2]*"FC"))+ylab(expression("-log"[10]*"FDR"))+
  scale_color_manual(values=c("dodgerblue3","gray50","firebrick3"))+
  guides(color=guide_legend(override.aes = list(size=1.5)))
p1

p2<-p1+geom_label(data=top_genes,mapping=aes(x=log2FoldChange,y=-log(padj,10),
                                             label=sample_ID),size=2)
p2

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

# heatmap
Genes_anno<-as.data.frame(top_genes$Expression)
rownames(Genes_anno)=rownames(top_genes)
colnames(Genes_anno)="Expression"
#prepare the gene data
countData_gene<-countData
countData_gene$gene<-rownames(countData)
top_count_seq<-left_join(top_genes,countData_gene,by=c("sample_ID"="gene"))
rownames(top_count_seq)<-top_genes$sample_ID
top_count_seq<-dplyr::select(top_count_seq,10:571)
norTop_count_seq<-log2(top_count_seq+1)

pheatmap::pheatmap(norTop_count_seq,
                   annotation_row = Genes_anno,
                   annotation_col = TN_Ident,
                   show_colnames = F,
                   annotation_names_row = F,
                   annotation_names_col = F)

# draw with complexheatmap
# make col_annotation
col_ha<-HeatmapAnnotation(df=TN_Ident,
                          show_annotation_name = F,
                          col = list(TN_Ident=c('tumor'='#D0104C','normal'='#0B346E')))
row_ha<-rowAnnotation(df=Genes_anno,
                      show_annotation_name = F,
                      col = list(Expression=c('Up-regulated'='#00896C','Down-regulated'='#F7C242')))

# high light the genes from Cox analysis
h_gl<-c("FN1","CFD","CDR2","DPP4","MMRN1")
index<-which(rownames(norTop_count_seq) %in% h_gl)
labs<-rownames(norTop_count_seq)[index]
lab2 = rowAnnotation(foo = anno_mark(at = index,
                                     labels = labs,
                                     lines_gp = gpar()))

p3<-Heatmap(as.matrix(norTop_count_seq),
        top_annotation=col_ha,
        left_annotation=row_ha,
        right_annotation=lab2,
        show_column_names=F,
        show_row_names=F)
p3

Heatmap(as.matrix(norTop_count_seq),
        top_annotation=col_ha,
        left_annotation=row_ha,
        right_annotation=lab2,
        show_column_names=F,
        show_row_names=F)

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
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)

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

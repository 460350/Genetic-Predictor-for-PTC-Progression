library(Seurat)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)

# visualization of MAPK13, STMN1, HMGB2
Idents(fhs_thyrocyte)<-'annotation'

FeaturePlot(fhs_thyrocyte,
            features = c('MAPK13','STMN1','HMGB2'),
            slot = 'data',
            min.cutoff = 1)
DimPlot(fhs_thyrocyte,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)

# identify the cells with marker gene expression
EXPR = GetAssayData(object=fhs_thyrocyte,assay="SCT",slot="data")["MAPK13",]
MAPK13_df=data.frame(MAPK13_positive= EXPR > 0)
fhs_thyrocyte<-AddMetaData(fhs_thyrocyte,MAPK13_df,col.name = 'MAPK13_pos')

EXPR = GetAssayData(object=fhs_thyrocyte,assay="SCT",slot="data")["STMN1",]
MAPK13_df=data.frame(MAPK13_positive= EXPR > 0)
fhs_thyrocyte<-AddMetaData(fhs_thyrocyte,MAPK13_df,col.name = 'STMN1_pos')

EXPR = GetAssayData(object=fhs_thyrocyte,assay="SCT",slot="data")["HMGB2",]
MAPK13_df=data.frame(MAPK13_positive= EXPR > 0)
fhs_thyrocyte<-AddMetaData(fhs_thyrocyte,MAPK13_df,col.name = 'HMGB2_pos')

summary(fhs_thyrocyte[[c('MAPK13_pos','STMN1_pos','HMGB2_pos')]])
# MAPK13_pos      STMN1_pos       HMGB2_pos      
# Mode :logical   Mode :logical   Mode :logical  
# FALSE:16031     FALSE:17050     FALSE:14086    
# TRUE :6495      TRUE :5476      TRUE :8440   

nrow(fhs_thyrocyte@meta.data[fhs_thyrocyte@meta.data[,39] == 'TRUE' & 
                               fhs_thyrocyte@meta.data[,40] == 'TRUE' & 
                               fhs_thyrocyte@meta.data[,41] == 'TRUE',]) # count rows with triple TRUE

# analyze the cells with triple positive
fhs_thyrocyte@meta.data$GL3T<-ifelse(fhs_thyrocyte@meta.data$MAPK13_pos == 'TRUE' & 
                                       fhs_thyrocyte@meta.data$STMN1_pos == 'TRUE' & 
                                       fhs_thyrocyte@meta.data$HMGB2_pos == 'TRUE', 2,
                                                ifelse(fhs_thyrocyte@meta.data$MAPK13_pos == 'FALSE' & 
                                                         fhs_thyrocyte@meta.data$STMN1_pos == 'FALSE' & 
                                                         fhs_thyrocyte@meta.data$HMGB2_pos == 'FALSE', 1,
                                                       0))
summary(as.factor(fhs_thyrocyte@meta.data$GL3T))
FeaturePlot(fhs_thyrocyte,
            features = 'GL3T',
            slot = 'data')

####GSEA Analysis####
fhs_thyrocyte@meta.data$GL3T[fhs_thyrocyte@meta.data$GL3T==2]<-'Triple_pos'
fhs_thyrocyte@meta.data$GL3T[fhs_thyrocyte@meta.data$GL3T==1]<-'Triple_neg'
fhs_thyrocyte@meta.data$GL3T[fhs_thyrocyte@meta.data$GL3T==0]<-'undefined'

table(fhs_thyrocyte@meta.data$GL3T)
# Triple_neg Triple_pos  undefined 
#   8964       1258      12304 

# DEGs between triple_pos and triple_neg
Idents(fhs_thyrocyte)<-'GL3T'
GL_DEGs<-FindMarkers(fhs_thyrocyte, ident.1="Triple_pos",ident.2="Triple_neg",
                                min.pct = 0.1,logfc.threshold = 0.25)
GL_DEGs$SYMBOL<-rownames(GL_DEGs)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
GL_DEGs<-GL_DEGs%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(GL_DEGs)<-GL_DEGs$SYMBOL

### enrichment analysis
geneList<-GL_DEGs[,2]
names(geneList)<-as.character(GL_DEGs[,6])
geneList<-sort(geneList,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
# select the KEGG pathway only
emC2_KEGG<-emC2
emC2_KEGG@result<-emC2_KEGG@result[grep("KEGG_", emC2@result$Description),]
enrichplot::dotplot(emC2_KEGG,showCategory = 20,split=".sign")+facet_grid(~.sign)
emC2_KEGG_mod<-emC2_KEGG
rownames(emC2_KEGG@result)
emC2_KEGG_mod@result<-emC2_KEGG@result[-c(2,3,4,5),]
enrichplot::dotplot(emC2_KEGG_mod,showCategory = 5,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(NA,0.94)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
emC5BP<-GSEA(geneList,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
emC5CC<-GSEA(geneList,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'Triple_pos & Triple_neg GSEA Result.RData')

####Analysis for the marker-gene high and normal group####
## define the marker genes high expression cells
expr = GetAssayData(object=fhs_thyrocyte,assay="SCT",slot="data")["MAPK13",]
summary(expr)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.2336  0.6931  2.0794 
MAPK13_expr=data.frame(MAPK13_high= expr > 1)
fhs_thyrocyte<-AddMetaData(fhs_thyrocyte,MAPK13_expr,col.name = 'MAPK13_high')

expr = GetAssayData(object=fhs_thyrocyte,assay="SCT",slot="data")["STMN1",]
summary(expr)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.0000  0.2268  0.0000  3.3673 
STMN1_expr=data.frame(STMN1_high= expr > 1)
fhs_thyrocyte<-AddMetaData(fhs_thyrocyte,STMN1_expr,col.name = 'STMN1_high')

expr = GetAssayData(object=fhs_thyrocyte,assay="SCT",slot="data")["HMGB2",]
summary(expr)
#  Min.    1st Qu.  Median    Mean  3rd Qu.    Max. 
#  0.0000  0.0000   0.0000  0.3456  0.6931  3.5553 
HMGB2_expr=data.frame(HMGB2_high= expr > 1)
fhs_thyrocyte<-AddMetaData(fhs_thyrocyte,HMGB2_expr,col.name = 'HMGB2_high')

summary(fhs_thyrocyte[[c('MAPK13_high','STMN1_high','HMGB2_high')]])

fhs_thyrocyte@meta.data$MAPK13_high[fhs_thyrocyte@meta.data$MAPK13_high=='TRUE']<-'High'
fhs_thyrocyte@meta.data$MAPK13_high[fhs_thyrocyte@meta.data$MAPK13_pos == 'FALSE'] <- "Zero"
fhs_thyrocyte@meta.data$MAPK13_high[fhs_thyrocyte@meta.data$MAPK13_high=='FALSE']<-'undefined'
table(fhs_thyrocyte@meta.data$MAPK13_high)

fhs_thyrocyte@meta.data$STMN1_high[fhs_thyrocyte@meta.data$STMN1_high=='TRUE']<-'High'
fhs_thyrocyte@meta.data$STMN1_high[fhs_thyrocyte@meta.data$STMN1_pos == 'FALSE'] <- "Zero"
fhs_thyrocyte@meta.data$STMN1_high[fhs_thyrocyte@meta.data$STMN1_high=='FALSE']<-'undefined'
table(fhs_thyrocyte@meta.data$STMN1_high)

fhs_thyrocyte@meta.data$HMGB2_high[fhs_thyrocyte@meta.data$HMGB2_high=='TRUE']<-'High'
fhs_thyrocyte@meta.data$HMGB2_high[fhs_thyrocyte@meta.data$HMGB2_pos == 'FALSE'] <- "Zero"
fhs_thyrocyte@meta.data$HMGB2_high[fhs_thyrocyte@meta.data$HMGB2_high=='FALSE']<-'undefined'
table(fhs_thyrocyte@meta.data$HMGB2_high)

# analyze the cells with triple positive
fhs_thyrocyte@meta.data$GL3H<-ifelse(fhs_thyrocyte@meta.data$MAPK13_high == 'High' & 
                                       fhs_thyrocyte@meta.data$STMN1_high == 'High' & 
                                       fhs_thyrocyte@meta.data$HMGB2_high == 'High', 2,
                                     ifelse(fhs_thyrocyte@meta.data$MAPK13_high == 'Zero' & 
                                              fhs_thyrocyte@meta.data$STMN1_high == 'Zero' & 
                                              fhs_thyrocyte@meta.data$HMGB2_high == 'Zero', 1,
                                            0))
fhs_thyrocyte@meta.data$GL3H[fhs_thyrocyte@meta.data$GL3H==2]<-'Triple_high'
fhs_thyrocyte@meta.data$GL3H[fhs_thyrocyte@meta.data$GL3H==1]<-'Triple_zero'
fhs_thyrocyte@meta.data$GL3H[fhs_thyrocyte@meta.data$GL3H==0]<-'Intermediate'
summary(as.factor(fhs_thyrocyte@meta.data$GL3H))
#Intermediate  Triple_high  Triple_zero 
#   13498           64         8964  

# DEGs between triple_high and triple_zero
Idents(fhs_thyrocyte)<-'GL3H'
GL3H_DEGs<-FindMarkers(fhs_thyrocyte, ident.1="Triple_high",ident.2=NULL,
                         min.pct = 0.1,logfc.threshold = 0.25)
GL3H_DEGs$SYMBOL<-rownames(GL3H_DEGs)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
GL3H_DEGs<-GL3H_DEGs%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(GL3H_DEGs)<-GL3H_DEGs$SYMBOL

### enrichment analysis
geneList<-GL3H_DEGs[,2]
names(geneList)<-as.character(GL3H_DEGs[,6])
geneList<-sort(geneList,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
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
enrichplot::dotplot(emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'Triple_high & Rest GSEA Result.RData')

# DEGs between MAPK13_high and MAPK13_zero
Idents(fhs_thyrocyte)<-'MAPK13_high'
MAPK13_DEGs<-FindMarkers(fhs_thyrocyte, ident.1="High",ident.2="Zero",
                     min.pct = 0.1,logfc.threshold = 0.25)
MAPK13_DEGs$SYMBOL<-rownames(MAPK13_DEGs)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
MAPK13_DEGs<-MAPK13_DEGs%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(MAPK13_DEGs)<-MAPK13_DEGs$SYMBOL

### enrichment analysis
geneList<-MAPK13_DEGs[,2]
names(geneList)<-as.character(MAPK13_DEGs[,6])
geneList<-sort(geneList,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term enriched under specific pvalueCutoff

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
# select the KEGG pathway only
emC2_KEGG<-emC2
emC2_KEGG@result<-emC2_KEGG@result[grep("KEGG_", emC2@result$Description),]
enrichplot::dotplot(emC2_KEGG,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term enriched under specific pvalueCutoff

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
emC5BP<-GSEA(geneList,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
emC5CC<-GSEA(geneList,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term enriched under specific pvalueCutoff

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'MAPK13_high & MAPK13_zero GSEA Result.RData')

# DEGs between STMN1_high and STMN1_zero
Idents(fhs_thyrocyte)<-'STMN1_high'
STMN1_DEGs<-FindMarkers(fhs_thyrocyte, ident.1="High",ident.2="Zero",
                         min.pct = 0.1,logfc.threshold = 0.25)
STMN1_DEGs$SYMBOL<-rownames(STMN1_DEGs)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
STMN1_DEGs<-STMN1_DEGs%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(STMN1_DEGs)<-STMN1_DEGs$SYMBOL

### enrichment analysis
geneList<-STMN1_DEGs[,2]
names(geneList)<-as.character(STMN1_DEGs[,6])
geneList<-sort(geneList,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
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
enrichplot::dotplot(emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'STMN1_high & STMN1_zero GSEA Result.RData')


# DEGs between HMGB2_high and HMGB2_zero
Idents(fhs_thyrocyte)<-'HMGB2_high'
HMGB2_DEGs<-FindMarkers(fhs_thyrocyte, ident.1="High",ident.2="Zero",
                        min.pct = 0.1,logfc.threshold = 0.25)
HMGB2_DEGs$SYMBOL<-rownames(HMGB2_DEGs)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
HMGB2_DEGs<-HMGB2_DEGs%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(HMGB2_DEGs)<-HMGB2_DEGs$SYMBOL

### enrichment analysis
geneList<-HMGB2_DEGs[,2]
names(geneList)<-as.character(HMGB2_DEGs[,6])
geneList<-sort(geneList,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
emH<-GSEA(geneList,TERM2GENE=H_t2g)
enrichplot::dotplot(emH,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term enriched under specific pvalueCutoff

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
emC2<-GSEA(geneList,TERM2GENE=C2_t2g)
enrichplot::dotplot(emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
# select the KEGG pathway only
emC2_KEGG<-emC2
emC2_KEGG@result<-emC2_KEGG@result[grep("KEGG_", emC2@result$Description),]
enrichplot::dotplot(emC2_KEGG,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is from KEGG in emC2

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
emC5BP<-GSEA(geneList,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
emC5CC<-GSEA(geneList,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is enriched under specific pvalueCutoff

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
emC5MF<-GSEA(geneList,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

# save the GSEA result
save(geneList,emC2,emC2_KEGG,emH,emC5BP,emC5CC,emC5MF,file = 'HMGB2_high & HMGB2_zero GSEA Result.RData')

# analysis of the percentage between cell subgroup and 3-Gene marker expression
table(Idents(fhs_thyrocyte),fhs_thyrocyte$GL3T)
table(Idents(fhs_thyrocyte),fhs_thyrocyte$GL3H)
# visualization


#### GSVA for scRNAseq Data ####
library(GSVA)
library(pheatmap)
# get the expression matrix of fhs_thyrocyte
expr<-as.data.frame(fhs_thyrocyte@assays$SCT@data)
expr<-as.matrix(expr)
meta<-fhs_thyrocyte@meta.data[,c('thca_ident','SCEVANpred')]
meta$thca_ident<-as.character(meta$thca_ident)
meta<-meta[,-2,drop=F]
# get the gene sets for gsva
m_df<-msigdbr(species = "Homo sapiens",category = "C2",subcategory = "CP:KEGG")
msigdbr_list<-split(x=m_df$gene_symbol,f=m_df$gs_name)
# run GSVA
kegg<-gsva(expr,msigdbr_list,kcdf="Gaussian",method='gsva',verbose=T)

# primary visualization with heatmap
pheatmap::pheatmap(kegg,show_rownames = T,show_colnames = F,annotation_col = meta,
                   fontsize_row = 5,filename = 'scRNA_gsva_heatmap.png',
                   width = 15,height = 12)

## cell type and pathway association
meta_cp<-meta%>%arrange(meta$thca_ident)
data<-kegg[,rownames(meta_cp)]
# get the summary of data
row_means<-apply(data,1,mean)
row_medians<-apply(data,1,median)

group<-factor(meta_cp[,'thca_ident'],ordered=F)
# visualization with heatmap
data1<-NULL

for(i in 0:(length(unique(group))-1)){
  ind<-which(group==i)
  dat<-apply(data[,ind],1,mean)
  data1<-cbind(data1,dat)
}
colnames(data1)<-c('HLA-','NMU-','NMB-','MTs-','MSMP-','IGFBP5-','RGS5-')
result<-t(scale(t(data1)))
result[result>2]=2
result[result<-2]=-2

pheatmap::pheatmap(result[1:20,],
                   cluster_rows=F,
                   cluster_cols=F,
                   show_rownames=T,
                   show_colnames=T,
                   color=colorRampPalette(c('blue','white','red'))(100),
                   cellwidth = 10,cellheight = 15,
                   fontsize = 10,
                   filename = 'gsva_celltype.png')

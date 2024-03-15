# The GSEA conducted here is totally based on unfiltered DEGs and new rank algorithm. The new rank algorithm elevate the power of p value in ranking.
# The previous one ranks was mostly based on fold change since the gene list were generated with cutoff pvalue <0.05 and relatively high percentage in cell clusters.
# Two figures are about GSEA results from two Scissor analyses and MAPK13, STMN1, HMGB2 three-gene model GSEA

library(tidyverse)
library(Seurat)
library(msigdbr)
library(clustree)
library(clusterProfiler)

# load the scRNAseq file
load('~/R Data/PTC_Pair/PTC without Normal Cells after CNV/combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData')

####Scissor analysis####
###stage###
Idents(fhs_thyrocyte)<-"Scissor_stage"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="1",ident.2="2",
                                min.pct = 0,logfc.threshold = 0)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL

# Rank score calculation
Scissor.de.markers$rankings<-sign(Scissor.de.markers$avg_log2FC)*(-log10(Scissor.de.markers$p_val))
write.csv(Scissor.de.markers,'Scissor_stagev2_ptc_de_genes_newgsea.csv')

# GSEA
### enrichment analysis
geneListT<-Scissor.de.markers[,7]
names(geneListT)<-as.character(Scissor.de.markers[,6])
geneListT<-sort(geneListT,decreasing = T)
plot(geneListT)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH<-GSEA(geneListT,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2<-GSEA(geneListT,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2,showCategory = 10,split=".sign")+facet_grid(~.sign)
t.emC2_KEGG<-t.emC2
t.emC2_KEGG@result<-t.emC2_KEGG@result[grep("KEGG_", t.emC2@result$Description),]
# no KEGG pathways are found in the C2 results
# try REACTOME pathways instead
t.emC2_REACTOME<-t.emC2
t.emC2_REACTOME@result<-t.emC2_REACTOME@result[grep("REACTOME_", t.emC2@result$Description),]
enrichplot::dotplot(t.emC2_REACTOME,showCategory = 9,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.15,0.75)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.15,0.55)

###DFS###
Idents(fhs_thyrocyte)<-"Scissor_DFS"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="1",ident.2="2",
                                min.pct = 0,logfc.threshold = 0)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL

# Rank score calculation
Scissor.de.markers$rankings<-sign(Scissor.de.markers$avg_log2FC)*(-log10(Scissor.de.markers$p_val))
write.csv(Scissor.de.markers,'Scissor_dfs_ptc_de_genes_newgsea.csv')
# GSEA
### enrichment analysis
geneListT<-Scissor.de.markers[,7]
names(geneListT)<-as.character(Scissor.de.markers[,6])
geneListT<-sort(geneListT,decreasing = T)
plot(geneListT)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH<-GSEA(geneListT,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH,showCategory = 5,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.15,0.63)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2<-GSEA(geneListT,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2,showCategory = 10,split=".sign")+facet_grid(~.sign)
t.emC2_KEGG<-t.emC2
t.emC2_KEGG@result<-t.emC2_KEGG@result[grep("KEGG_", t.emC2@result$Description),]
enrichplot::dotplot(t.emC2_KEGG,showCategory = 6,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.17,0.55)
# try REACTOME pathways instead
t.emC2_REACTOME<-t.emC2
t.emC2_REACTOME@result<-t.emC2_REACTOME@result[grep("REACTOME_", t.emC2@result$Description),]
enrichplot::dotplot(t.emC2_REACTOME,showCategory = 6,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.1,0.55)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

####Three Gene Model####
Idents(fhs_thyrocyte)<-"GL3T"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="Triple_pos",ident.2="Triple_neg",
                                min.pct = 0,logfc.threshold = 0)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL

# Rank score calculation
Scissor.de.markers$rankings<-sign(Scissor.de.markers$avg_log2FC)*(-log10(Scissor.de.markers$p_val))
write.csv(Scissor.de.markers,'Scissor_3Gene_ptc_de_genes_newgsea.csv')
# GSEA
### enrichment analysis
geneListT<-Scissor.de.markers[,7]
names(geneListT)<-as.character(Scissor.de.markers[,6])
geneListT<-sort(geneListT,decreasing = T)
plot(geneListT)
# The rankings of STMN1, PASK, PPP1R7, HMGB2, MAPK13, XIST, NMU, KLK11 are 'Inf' due to extremely small p value.
# Assign a relatively small p value to it. The p value is arbitarily set as 10^(-300) since they are indeed very different
# and it is the smallest number in r in a data.frame
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'STMN1'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'PASK'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'PPP1R7'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'HMGB2'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'MAPK13'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'XIST'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'NMU'] <- 10^(-300)
Scissor.de.markers$p_val[rownames(Scissor.de.markers) == 'KLK11'] <- 10^(-300)
# Rank score calculation again
Scissor.de.markers$rankings<-sign(Scissor.de.markers$avg_log2FC)*(-log10(Scissor.de.markers$p_val))
write.csv(Scissor.de.markers,'Scissor_3Gene_ptc_de_genes_newgsea.csv')

geneListT<-Scissor.de.markers[,7]
names(geneListT)<-as.character(Scissor.de.markers[,6])
geneListT<-sort(geneListT,decreasing = T)
plot(geneListT)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH<-GSEA(geneListT,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH,showCategory = 5,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2<-GSEA(geneListT,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2,showCategory = 10,split=".sign")+facet_grid(~.sign)
t.emC2_KEGG<-t.emC2
t.emC2_KEGG@result<-t.emC2_KEGG@result[grep("KEGG_", t.emC2@result$Description),]
enrichplot::dotplot(t.emC2_KEGG,showCategory = 5,split=".sign")+facet_grid(~.sign)
# try REACTOME pathways instead
t.emC2_REACTOME<-t.emC2
t.emC2_REACTOME@result<-t.emC2_REACTOME@result[grep("REACTOME_", t.emC2@result$Description),]
enrichplot::dotplot(t.emC2_REACTOME,showCategory = 5,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

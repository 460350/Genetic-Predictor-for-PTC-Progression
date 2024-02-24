library(tidyverse)
library(Seurat)
library(harmony)
library(clustree)
library(msigdbr)
library(enrichplot)
library(clusterProfiler)

####merge####
# load the files
fs_ptcTN<-merge(x = ptc2_3tn,
                y = ptc5t,
                project = "thyroidT")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc10t,
                project = "thyroidT")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc5rin,
                project = "thyroidN")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc6rin,
                project = "thyroidN")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc7rin,
                project = "thyroidN")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc10rin,
                project = "thyroidN")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc11rin,
                project = "thyroidN")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc4sc,
                project = "thyroidSC")
fs_ptcTN<-merge(x = fs_ptcTN,
                y = ptc11sc,
                project = "thyroidSC")
####
####normalization and feature selection####
# findvariablefeatures has been done before, so skip this step
fs_ptcTN<-fs_ptcTN%>%
  SCTransform(vars.to.regress = c("mitoRatio","riboRatio","CC.diff"))
ElbowPlot(fs_ptcTN, ndims = 50)
# calculate PCs using variable features determined by SCTransform
fs_ptcTN<-RunPCA(fs_ptcTN,assay = "SCT",npcs = 50)
ElbowPlot(fs_ptcTN, ndims = 50)
####

####integration####
fhs_ptcTN<-RunHarmony(fs_ptcTN,
                      group.by.vars = "seq_folder",
                      reduction = "pca",assay.use = "SCT",reduction.save = "harmony")
####
####clustering & visualization####
fhs_ptcTN<-RunUMAP(object = fhs_ptcTN,reduction = "harmony",assay = "SCT",dims = 1:40)%>%
  FindNeighbors(reduction = "harmony")%>%
  FindClusters(resolution = c(0.2,0.4,0.6,0.8,1.0))

clustree::clustree(fhs_ptcTN@meta.data,prefix="SCT_snn_res.")

Idents(object = fhs_ptcTN)<-"SCT_snn_res.0.4"

fhs_ptcTN<-RunUMAP(object = fhs_ptcTN,reduction = "harmony",assay = "SCT",dims = 1:40)
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        group.by = "seq_folder",
        split.by = "origin")
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        split.by = "seq_folder",
        ncol = 2)
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        group.by = "SCEVANpred")
####
####Overall QC####
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio","riboRatio")
FeaturePlot(fhs_ptcTN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)
RidgePlot(fhs_ptcTN, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
####
####SAVE####
save(fhs_ptcTN,file="PTC Sample/combinedCNV_ptc_TNSC.RData")
####
####annotation####
annotations<-read.csv("annotation.csv")

fhsptcTN_markers<-FindAllMarkers(object = fhs_ptcTN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Combine marker with gene description
fhsptcTN_markers<-fhsptcTN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 20 markers per cluster
fhsptcTN_top20<-fhsptcTN_markers%>%
  group_by(cluster)%>%
  top_n(n = 20, wt = avg_log2FC)

# Plot the specific marker distribution in cluster
FeaturePlot(object = fhs_ptcTN, 
            features = c("TG"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE,
            split.by = "origin")
VlnPlot(object = fhs_ptcTN,
        features = c("TG","TPO"))

####original cluster 9 needs to be removed due to high immune marker expression####
VlnPlot(fhs_ptcTN,
        features = c('IGKC','PTPRC','TRBC2','LYZ'),
        ncol = 2)
fhs_ptcTN<-subset(fhs_ptcTN,idents = "9",invert=T)
####


#####

####reset the cluster identity to default####
Idents(object = fhs_ptcTN)<-"SCT_snn_res.0.4"


####cluster identity for 0.4####
fhs_ptcTN<-RenameIdents(object = fhs_ptcTN,
                        "0"="NMB-Thyrocyte",
                        "1"="HLA-Thyrocyte",
                        "2"="NMU-Thyrocyte",
                        "3"="CAF",
                        "4"="Endothelial cell",
                        "5"="Fibroblast",
                        "6"="IGFBP5-Thyrocyte",
                        "7"="MSMP-Thyrocyte",
                        "8"="MTs-Thyrocyte",
                        "9"="IGFBP5-Thyrocyte",
                        "10"="RGS5-Thyrocyte",
                        "11"="NMB-Thyrocyte")
p1<-DimPlot(fhs_ptcTN,
            reduction = "umap",
            label = T,
            label.size = 4,
            repel = T)
p1

DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        split.by = "seq_folder")
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 3,
        repel = T,
        split.by = "origin")
# SCEVANpred result
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        group.by = "SCEVANpred")
# evaluate the result with several tumor markers
FeaturePlot(object = fhs_ptcTN,
            features = c("TG","EPCAM"),
            order = T,
            min.cutoff = 'q10',
            label = T,
            repel = T,
            split.by = "origin")
# evaluate the contribution of each sample
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        split.by = "seq_folder",
        ncol = 4)
DimPlot(fhs_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        split.by = "seq_folder",
        group.by = "SCEVANpred",
        ncol = 4)

ptc_features<-c("TSHR","TG","TPO","EPCAM","KRT19","NMB","CXCL14","HLA-DRA","NMU","IGFBP5",
                "MSMP","MT1X","CXCL2","RGS5","ACTA2","COL1A2","COL1A1","VWF","PECAM1")

p2<-DotPlot(fhs_ptcTN,features = ptc_features,cols = c("blue","red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
# change the order of x-axis in p2
nonImmune.order.levels<-c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
                          "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte","CAF",
                          "Fibroblast","Endothelial cell")
Idents(fhs_ptcTN)<-factor(Idents(fhs_ptcTN),levels = nonImmune.order.levels)
p2<-DotPlot(fhs_ptcTN,features = ptc_features,cols = c("blue","red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
p2
#####

####percentage of each cluster among the thyrocytes and subgroup comparison####
fhs_thyrocyte<-subset(fhs_ptcTN,
                      idents=c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
                               "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte"))
## SAVE
save(fhs_thyrocyte,file="PTC Sample/combinedCNV_ptc_TNSC_thyrocyte.RData")
# check the contribution of each sample
table(fhs_thyrocyte$orig.ident)
# percentage of each identified cluster
prop.table(table(Idents(fhs_thyrocyte)))
# make a table based on the results above
table(Idents(fhs_thyrocyte),fhs_thyrocyte$orig.ident)
prop.table(table(Idents(fhs_thyrocyte),fhs_thyrocyte$orig.ident),margin = 2)

# calculation
cellRatio<-prop.table(table(Idents(fhs_thyrocyte),fhs_thyrocyte$orig.ident),margin = 2)%>%
  as.data.frame()
cellRatio
colorCount<-length(unique(cellRatio$Var1))
ggplot(cellRatio)+
  geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour='#222222')+
  theme_classic()+
  labs(x='Sample',y='Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))

# calculation
CellRatio<-prop.table(table(Idents(fhs_thyrocyte),fhs_thyrocyte$origin),margin = 2)%>%
  as.data.frame()
CellRatio
colorCount<-length(unique(CellRatio$Var1))
ggplot(CellRatio)+
  geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour='#222222')+
  theme_classic()+
  labs(x='Sample',y='Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))

## make comparison between tumor and lymph node quantitively
library(reshape2)
cellper<-dcast(cellRatio%>%as_tibble()%>%filter(.,! Var2 %in% "sc")%>%as.data.frame(),
               Var2~Var1,value.var = "Freq")
rownames(cellper)<-cellper[,1]
cellper<-cellper[,-1]

sample<-c('ptc2t','ptc2len','ptc3t','ptc3len','ptc5t','ptc5rin','ptc10t','ptc10rin','ptc6rin','ptc7rin','ptc11rin',
          'ptc4sc','ptc11sc')
group<-c('tumor','lymph node','tumor','lymph node','tumor','lymph node','tumor','lymph node',
         'lymph node','lymph node','lymph node','subcutaneous','subcutaneous')
# create the data.frame
samples<-data.frame(sample,group)
rownames(samples)<-samples$sample

cellper$sample<-samples[rownames(cellper),'sample']
cellper$group<-samples[rownames(cellper),'group']

# figure production
pplist=list()
thca_group<-c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
              "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte")
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in thca_group){
  cellper_ = cellper %>% dplyr::select(one_of(c('sample','group',group_)))# choose one line of data
  colnames(cellper_)=c('sample','group','percent')# name the chosen data
  cellper_$percent = as.numeric(cellper_$percent)# numeric data
  cellper_<-cellper_ %>% group_by(group) %>% mutate(upper = quantile(percent,0.75),
                                                    lower = quantile(percent,0.25),
                                                    mean = mean(percent),
                                                    median = median(percent))
  print(group_)
  print(cellper_$median)
  # plotting
  pp1 = ggplot(cellper_,aes(x=group,y=percent))+ 
    geom_jitter(shape=21,aes(fill=group),width = 0.25)+
    stat_summary(fun=mean,geom = "point",color='grey60')+
    theme_cowplot()+
    theme(axis.text = element_text(size=10),axis.text.x = element_text(angle = 45, hjust = 1),axis.title = element_text(size=10),
          legend.text = element_text(size=10),legend.title = element_text(size=10),
          plot.title = element_text(size=10,face='plain'),legend.position = 'none')+
    labs(title = group_,y='Percentage')+
    geom_errorbar(aes(ymin = lower, ymax = upper),col='grey60',width=1)
  # t test
  labely = max(cellper_$percent)
  compare_means(percent~group,data=cellper_)
  my_comparisons<-list(c("lymph node","subcutaneous"))
  pp1 = pp1+stat_compare_means(comparisons=my_comparisons,size=2,method='t.test')
  pplist[[group_]] = pp1
}

library(cowplot)
plot_grid(pplist[['NMB-Thyrocyte']],
          pplist[['HLA-Thyrocyte']],
          pplist[['NMU-Thyrocyte']],
          pplist[['IGFBP5-Thyrocyte']],
          pplist[['MSMP-Thyrocyte']],
          pplist[['MTs-Thyrocyte']],
          pplist[['RGS5-Thyrocyte']])
#####

####Hypoxia Score####
H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
# prepare the emt genes
emt<-H_t2g%>%
  filter(gs_name=="HALLMARK_HYPOXIA")
emt<-as.data.frame(emt)
emt<-emt$gene_symbol
# check the overlap between emt gene list and current gene list from sample
emt<-intersect(emt,rownames(GetAssayData(fhs_ptcTN,slot = "scale.data")));length(emt)
emt_features<-list(emt)
fhs_ptcTN <- AddModuleScore(
  object = fhs_ptcTN,
  features = emt_features,
  ctrl = 100,
  name = 'Hypoxia_score'
)

# significant
library(ggsignif)
fhs_ptcTN@meta.data$annotation<-Idents(fhs_ptcTN)
p4 <- ggplot(subset(x=fhs_ptcTN, 
                    idents=c("IGFBP5-Thyrocyte","MTs-Thyrocyte","NMB-Thyrocyte","NMU-Thyrocyte",
                             "HLA-Thyrocyte","MSMP-Thyrocyte","RGS5-Thyrocyte"))@meta.data, 
             aes(x=annotation, y=Hypoxia_score1, color=annotation)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("IGFBP5-Thyrocyte","MSMP-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","MTs-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","RGS5-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","NMU-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","HLA-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","NMB-Thyrocyte")),
              y_position = c(0.56,0.60,0.64,0.68,0.72,0.76,0.78,0.82),
              map_signif_level = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none") + xlab(NULL) + ylab("Hypoxia Scores")
p4

VlnPlot(fhs_ptcTN,features = "Hypoxia_score1")
VlnPlot(subset(x=fhs_ptcTN, 
               idents=c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
                        "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte")),
        features = "Hypoxia_score1")

FeaturePlot(object = fhs_ptcTN, 
            features = c("Hypoxia_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)

FeaturePlot(object = subset(x=fhs_ptcTN, 
                            idents=c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
                                     "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte")), 
            features = c("EMT_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)
#####

####TDS calculation####
TDS_features<-c("NKX2-1","DUOX1","DUOX2","PAX8","SCL26A4","FOXE1","TG","TSHR",
                "THRA","DIO2","GLIS3","TPO","DIO1","THRB","SLC5A5","SLC5A8")
TDS_features<-list(TDS_features)

fhs_ptcTN <- AddModuleScore(
  object = fhs_ptcTN,
  features = TDS_features,
  ctrl = 100,
  name = 'TDS_score'
)

# significant
library(ggsignif)
fhs_ptcTN@meta.data$annotation<-Idents(fhs_ptcTN)
p5 <- ggplot(subset(x=fhs_ptcTN, 
                    idents=c("IGFBP5-Thyrocyte","MTs-Thyrocyte","NMB-Thyrocyte","NMU-Thyrocyte",
                             "HLA-Thyrocyte","MSMP-Thyrocyte","RGS5-Thyrocyte"))@meta.data, 
             aes(x=annotation, y=TDS_score1, color=annotation)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("IGFBP5-Thyrocyte","MSMP-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","MTs-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","RGS5-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","NMU-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","HLA-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","NMB-Thyrocyte")),
              y_position = c(0.94,1.00,1.06,1.12,1.18,1.24,1.30,1.36),
              map_signif_level = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + xlab(NULL) + ylab("TDS Scores")
p5

VlnPlot(fhs_ptcTN,features = "TDS_score1")
VlnPlot(subset(x=fhs_ptcTN, 
               idents=c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
                        "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte")),
        features = "TDS_score1")

FeaturePlot(object = fhs_ptcTN, 
            features = c("TDS_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)

FeaturePlot(object = subset(x=fhs_ptcTN, 
                            idents=c("NMB-Thyrocyte","HLA-Thyrocyte","NMU-Thyrocyte","IGFBP5-Thyrocyte",
                                     "MSMP-Thyrocyte","MTs-Thyrocyte","RGS5-Thyrocyte")), 
            features = c("TDS_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)
#####
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
# Lineage comparison & Call DEGs
## IGFBP5 & RGS5
# pre-filter features that are detected at <50% frequency in both group, 
# that have less than tow-old change between the average expression, 
# that detection percentage across the two groups are similar (within 0.25)
Tumor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="RGS5-Thyrocyte",
                              min.pct = 0.5,logfc.threshold = log(2),min.diff.pct = 0.25)
head(Tumor.de.markers)
Tumor.de.markers$SYMBOL<-rownames(Tumor.de.markers)
# 0 outcome

Tumor.de.markers.more<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="RGS5-Thyrocyte")
head(Tumor.de.markers.more)
Tumor.de.markers.more$SYMBOL<-rownames(Tumor.de.markers.more)

ggplot(data = Tumor.de.markers.more,aes(x=-log10(p_val),y=-log10(p_val_adj)))+
  geom_point()+
  geom_smooth(method = lm)
# further remove ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Tumor.de.markers.more<-Tumor.de.markers.more%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Tumor.de.markers.more)<-Tumor.de.markers.more$SYMBOL
Tumor.de.markers.more<-Tumor.de.markers.more%>%
  dplyr::filter(p_val_adj<0.05)
# export the gene list
write.csv(Tumor.de.markers.more,'IGFBP5_RGS5_ptc_de_genes.csv')

### enrichment analysis
geneListT.m<-Tumor.de.markers.more[,2]
names(geneListT.m)<-as.character(Tumor.de.markers.more[,6])
geneListT.m<-sort(geneListT.m,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH.m<-GSEA(geneListT.m,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH.m,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2.m<-GSEA(geneListT.m,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2.m,showCategory = 20,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT.m,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is selected

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT.m,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT.m,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is selected
## KEGG with online method
## output gene list
write.table(names(geneListT.m),'exportedData/IGFBP5_RGS5_geneListTm.csv',sep = "",row.names = F)
## data analyzed at http://kobas.cbi.pku.edu.cn/

## IGFBP5 & MSMP (the KEGG analysis is good)
# pre-filter features that are detected at <50% frequency in both group, 
# that have less than tow-old change between the average expression, 
# that detection percentage across the two groups are similar (within 0.25)
Tumor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="MSMP-Thyrocyte",
                              min.pct = 0.5,logfc.threshold = log(2),min.diff.pct = 0.25)
head(Tumor.de.markers)
Tumor.de.markers$SYMBOL<-rownames(Tumor.de.markers)

Tumor.de.markers.more<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="MSMP-Thyrocyte")
head(Tumor.de.markers.more)
Tumor.de.markers.more$SYMBOL<-rownames(Tumor.de.markers.more)
# further remove ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Tumor.de.markers.more<-Tumor.de.markers.more%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Tumor.de.markers.more)<-Tumor.de.markers.more$SYMBOL
Tumor.de.markers.more<-Tumor.de.markers.more%>%
  dplyr::filter(p_val_adj<0.05)
write.csv(Tumor.de.markers.more,'IGFBP5_MSMP_ptc_de_genes.csv')
### enrichment analysis
geneListT<-Tumor.de.markers[,2]
names(geneListT)<-as.character(Tumor.de.markers[,6])
geneListT<-sort(geneListT,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH<-GSEA(geneListT,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is selected

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2<-GSEA(geneListT,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is selected

### more
geneListT.m<-Tumor.de.markers.more[,2]
names(geneListT.m)<-as.character(Tumor.de.markers.more[,6])
geneListT.m<-sort(geneListT.m,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH.m<-GSEA(geneListT.m,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH.m,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is selected

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2.m<-GSEA(geneListT.m,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2.m,showCategory = 20,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT.m,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT.m,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)
# no term is selected

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT.m,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is selected

## KEGG with online method
## output gene list
write.table(names(geneListT),'exportedData/IGFBP5_MSMP_geneListT.csv',sep = "",row.names = F)
write.table(names(geneListT.m),'exportedData/IGFBP5_MSMP_geneListTm.csv',sep = "",row.names = F)
## data analyzed at http://kobas.cbi.pku.edu.cn/

## IGFBP5 & HLA (the KEGG analysis is good)
# pre-filter features that are detected at <50% frequency in both group, 
# that have less than tow-old change between the average expression, 
# that detection percentage across the two groups are similar (within 0.25)
Tumor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="HLA-Thyrocyte",
                              min.pct = 0.5,logfc.threshold = log(2),min.diff.pct = 0.25)
head(Tumor.de.markers)
Tumor.de.markers$SYMBOL<-rownames(Tumor.de.markers)

Tumor.de.markers.more<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="HLA-Thyrocyte")
head(Tumor.de.markers.more)
Tumor.de.markers.more$SYMBOL<-rownames(Tumor.de.markers.more)
# further remove ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Tumor.de.markers.more<-Tumor.de.markers.more%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Tumor.de.markers.more)<-Tumor.de.markers.more$SYMBOL
Tumor.de.markers.more<-Tumor.de.markers.more%>%
  dplyr::filter(p_val_adj<0.05)
write.csv(Tumor.de.markers.more,'IGFBP5_HLA_ptc_de_genes.csv')
### enrichment analysis
geneListT<-Tumor.de.markers[,2]
names(geneListT)<-as.character(Tumor.de.markers[,6])
geneListT<-sort(geneListT,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH<-GSEA(geneListT,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is selected

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2<-GSEA(geneListT,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is selected

### more
geneListT.m<-Tumor.de.markers.more[,2]
names(geneListT.m)<-as.character(Tumor.de.markers.more[,6])
geneListT.m<-sort(geneListT.m,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH.m<-GSEA(geneListT.m,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH.m,showCategory = 20,split=".sign")+facet_grid(~.sign)

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2.m<-GSEA(geneListT.m,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2.m,showCategory = 20,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT.m,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT.m,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT.m,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is selected
## KEGG with online method
## output gene list
write.table(names(geneListT),'exportedData/IGFBP5_HLA_geneListT.csv',sep = "",row.names = F)
write.table(names(geneListT.m),'exportedData/IGFBP5_HLA_geneListTm.csv',sep = "",row.names = F)
## data analyzed at http://kobas.cbi.pku.edu.cn/
#####

## IGFBP5 & NMU
# pre-filter features that are detected at <50% frequency in both group, 
# that have less than tow-old change between the average expression, 
# that detection percentage across the two groups are similar (within 0.25)
Tumor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="NMU-Thyrocyte",
                              min.pct = 0.5,logfc.threshold = log(2),min.diff.pct = 0.25)
head(Tumor.de.markers)
Tumor.de.markers$SYMBOL<-rownames(Tumor.de.markers)
# 14 outcome

Tumor.de.markers.more<-FindMarkers(fhs_thyrocyte, ident.1="IGFBP5-Thyrocyte",ident.2="NMU-Thyrocyte")
head(Tumor.de.markers.more)
Tumor.de.markers.more$SYMBOL<-rownames(Tumor.de.markers.more)
# 430 outcome

ggplot(data = Tumor.de.markers.more,aes(x=-log10(p_val),y=-log10(p_val_adj)))+
  geom_point()+
  geom_smooth(method = lm)
# further remove ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Tumor.de.markers.more<-Tumor.de.markers.more%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Tumor.de.markers.more)<-Tumor.de.markers.more$SYMBOL
Tumor.de.markers.more<-Tumor.de.markers.more%>%
  dplyr::filter(p_val_adj<0.05)
# 365 outcome
# export the gene list
write.csv(Tumor.de.markers.more,'IGFBP5_NMU_ptc_de_genes.csv')

### enrichment analysis
geneListT.m<-Tumor.de.markers.more[,2]
names(geneListT.m)<-as.character(Tumor.de.markers.more[,6])
geneListT.m<-sort(geneListT.m,decreasing = T)

H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
t.emH.m<-GSEA(geneListT.m,TERM2GENE=H_t2g)
enrichplot::dotplot(t.emH.m,showCategory = 20,split=".sign")+facet_grid(~.sign)
# no term is selected

C2_t2g<-msigdbr(species = "Homo sapiens",category="C2")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC2.m<-GSEA(geneListT.m,TERM2GENE=C2_t2g)
enrichplot::dotplot(t.emC2.m,showCategory = 20,split=".sign")+facet_grid(~.sign)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT.m,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT.m,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT.m,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)

## KEGG with online method
## output gene list
write.table(names(geneListT),'exportedData/IGFBP5_NMU_geneListT.csv',sep = "",row.names = F)
write.table(names(geneListT.m),'exportedData/IGFBP5_NMU_geneListTm.csv',sep = "",row.names = F)
## data analyzed at http://kobas.cbi.pku.edu.cn/
library(tidyverse)
library(Seurat)
library(msigdbr)
library(clustree)
library(clusterProfiler)

# load the non-immune cells sample
load('combinedCNV_ptc_TNSC_thyrocyte.RData')

# remove the cells identified as normal cells
Idents(fhs_thyrocyte)<-"SCEVANpred"
fhs_thyrocyte<-subset(fhs_thyrocyte,
                      idents=c("tumor","filtered"))

# check the result after normal cells removal
Idents(fhs_thyrocyte)<-"annotation"
DimPlot(fhs_thyrocyte)

# new dimensional reduction after subsetting the original group
fhs_thyrocyte<-RunUMAP(object = fhs_thyrocyte,reduction = "harmony",assay = "SCT",dims = 1:40)%>%
  FindNeighbors(reduction = "harmony")%>%
  FindClusters(resolution = c(0.2,0.4,0.6,0.8,1.0))

clustree::clustree(fhs_thyrocyte@meta.data,prefix="SCT_snn_res.")

Idents(object = fhs_thyrocyte)<-"SCT_snn_res.0.2"

fhs_thyrocyte<-RunUMAP(object = fhs_thyrocyte,reduction = "harmony",assay = "SCT",dims = 1:40)
DimPlot(fhs_thyrocyte,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)
DimPlot(fhs_thyrocyte,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        group.by = 'SCT_snn_res.0.2')
# the difference between SCT_snn_res. and original annotation is that it further split the IGFBP5 group, separating
# a minor subgroup that is different from major one (probably more benign). Keep the original Idents
Idents(fhs_thyrocyte)<-"annotation"

# SAVE the result without normal cells
save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

# the percentage of each sample
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
emt<-intersect(emt,rownames(GetAssayData(fhs_thyrocyte,slot = "scale.data")));length(emt)
emt_features<-list(emt)
fhs_thyrocyte <- AddModuleScore(
  object = fhs_thyrocyte,
  features = emt_features,
  ctrl = 100,
  name = 'Hypoxia_score'
)

# significant
library(ggsignif)
fhs_thyrocyte@meta.data$annotation<-Idents(fhs_thyrocyte)
p1 <- ggplot(fhs_thyrocyte@meta.data, 
             aes(x=annotation, y=Hypoxia_score1, color=annotation)) + 
  geom_boxplot() +
  geom_signif(comparisons = list(c("IGFBP5-Thyrocyte","MSMP-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","MTs-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","RGS5-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","NMU-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","HLA-Thyrocyte"),
                                 c("IGFBP5-Thyrocyte","NMB-Thyrocyte")),
              y_position = c(0.60,0.64,0.68,0.72,0.76,0.80,0.84,0.88),
              map_signif_level = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + xlab(NULL) + ylab("Hypoxia Scores")
p1

VlnPlot(fhs_thyrocyte,features = "Hypoxia_score1")

FeaturePlot(object = fhs_thyrocyte, 
            features = c("Hypoxia_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)

# calculate the exact median value of each subgroup
fhs_thyrocyte@meta.data %>%
  group_by(annotation) %>%
  summarise(median_pts = median(Hypoxia_score1, na.rm=TRUE))
#####

####TDS calculation####
TDS_features<-c("NKX2-1","DUOX1","DUOX2","PAX8","SCL26A4","FOXE1","TG","TSHR",
                "THRA","DIO2","GLIS3","TPO","DIO1","THRB","SLC5A5","SLC5A8")
TDS_features<-list(TDS_features)

fhs_thyrocyte <- AddModuleScore(
  object = fhs_thyrocyte,
  features = TDS_features,
  ctrl = 100,
  name = 'TDS_score'
)

# significant
library(ggsignif)
fhs_thyrocyte@meta.data$annotation<-Idents(fhs_thyrocyte)
p2 <- ggplot(fhs_thyrocyte@meta.data, 
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
p2

VlnPlot(fhs_thyrocyte,features = "TDS_score1")

FeaturePlot(object = fhs_thyrocyte, 
            features = c("TDS_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)

# calculate the exact median value of each subgroup
fhs_thyrocyte@meta.data %>%
  group_by(annotation) %>%
  summarise(median_pts = median(TDS_score1, na.rm=TRUE))
#####

####EMT Score####
H_t2g<-msigdbr(species = "Homo sapiens",category="H")%>%
  dplyr::select(gs_name,gene_symbol)
# prepare the emt genes
emt<-H_t2g%>%
  filter(gs_name=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
emt<-as.data.frame(emt)
emt<-emt$gene_symbol
# check the overlap between emt gene list and current gene list from sample
emt<-intersect(emt,rownames(GetAssayData(fhs_thyrocyte,slot = "scale.data")));length(emt)
emt_features<-list(emt)
fhs_thyrocyte <- AddModuleScore(
  object = fhs_thyrocyte,
  features = emt_features,
  ctrl = 100,
  name = 'EMT_score'
)

# significant
library(ggsignif)
fhs_thyrocyte@meta.data$annotation<-Idents(fhs_thyrocyte)
p3 <- ggplot(fhs_thyrocyte@meta.data, 
             aes(x=annotation, y=EMT_score1, color=annotation)) + 
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
        legend.position = "none") + xlab(NULL) + ylab("EMT Scores")
p3

VlnPlot(fhs_thyrocyte,features = "EMT_score1")

FeaturePlot(object = fhs_thyrocyte, 
            features = c("EMT_score1"),
            order = TRUE, 
            label = TRUE,
            repel = TRUE)

# calculate the exact median value of each subgroup
fhs_thyrocyte@meta.data %>%
  group_by(annotation) %>%
  summarise(median_pts = median(EMT_score1, na.rm=TRUE))
#####

####slingshot analysis####
library(slingshot)
library(scater)
library(scran)
library(grDevices)
library(RColorBrewer)
library(tradeSeq)
library(ggbeeswarm)
library(ComplexHeatmap)
library(magick)

# load the seurat object
# save the object as separate matrices for input in slingshot
dimred<-fhs_thyrocyte@reductions$umap@cell.embeddings
fhs_thyrocyte$thca_ident<-Idents(fhs_thyrocyte)
clustering<-fhs_thyrocyte$thca_ident

# add ARPC1B and MYL12B in it since they are not among the highly-variable genes but significant 
# in predicting TCGA_THCA overall survival
m_counts<-as.matrix(fhs_thyrocyte@assays$RNA@counts[c(fhs_thyrocyte@assays$SCT@var.features,'ARPC1B','MYL12B'),])

# trajectory inference with slingshot
lineages<-getLineages(data = dimred,clusterLabels = clustering,
                      start.clus = "IGFBP5-Thyrocyte")
lineages

# Plot the lineages
# Define a color pallete to use
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(lineages), lwd = 3, col = "black")
opar<-par(no.readonly=TRUE)
par(opar)

# defining principal curves
curves <- getCurves(lineages,
                    stretch = 0,
                    extend = 'n',
                    allow.breaks = FALSE, shrink = T)
curves

curves@metadata$lineages

plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(SlingshotDataSet(curves), lwd = 3, col = "black")

# finding differentially expressed genes
library(tradeSeq)
# use highly variable genes to speed up the analysis
seu_variable_genes<-VariableFeatures(fhs_thyrocyte)[1:2000]

seu_variable_genes<-seu_variable_genes[!grepl('^RP[LS]', seu_variable_genes)]
seu_variable_genes<-seu_variable_genes[!grepl('^IG[HL]', seu_variable_genes)]
seu_variable_genes<-seu_variable_genes[!grepl('^MT-', seu_variable_genes)]

# add ARPC1B and MYL12B in it since they are not among the highly-variable genes but significant 
# in predicting TCGA_THCA overall survival
m_seu_sce<-fitGAM(counts = m_counts, sds = curves,
                  genes = c(seu_variable_genes,'ARPC1B','MYL12B'))
save(m_seu_sce,curves,m_counts,file="Seu_sling_trade_Thyrocyte_m.RData")

table(rowData(m_seu_sce)$tradeSeq$converged)

plotGeneCount(curves,m_counts,clusters=clustering,models=m_seu_sce)

# add slingshot result to seurat oject
fhs_thyrocyte$slingPseudotime_1<-m_seu_sce$crv$pseudotime.Lineage1
fhs_thyrocyte$slingPseudotime_2<-m_seu_sce$crv$pseudotime.Lineage2
fhs_thyrocyte$slingPseudotime_3<-m_seu_sce$crv$pseudotime.Lineage3

# SAVE the result without normal cells
save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

slingshot_df<-fhs_thyrocyte@meta.data
rm(fhs_thyrocyte)
# define a function to plot
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% dplyr::filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, m_counts, gene = feature_id[1], clusters = clustering, models = m_seu_sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(m_seu_sce, as.matrix(m_counts), gene = feature_id[1]))
}

# genes that change with pseudotime
pseudotime_association <- associationTest(m_seu_sce,lineages = T)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

cowplot::plot_grid(plotGeneCount(curves, m_counts, gene = "CD74", clusters = clustering, models = m_seu_sce) + ggplot2::theme(legend.position = "none"), 
                   plotSmoothers(m_seu_sce, as.matrix(m_counts), gene = "CD74"))

                    
# heatmap for gene expression during the pseudotime
heat_count<-fhs_thyrocyte@assays$RNA@counts
slingshot_df<-fhs_thyrocyte@meta.data
# lineage 1
ggplot(slingshot_df, aes(x = slingPseudotime_1, y = thca_ident,
                         colour = thca_ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime 1")

top_pa_genes1<-top_pa_genes
# Replace all NA values with 1 so those genes won't be excluded in the following filtering
top_pa_genes1[is.na(top_pa_genes1)] <- 1
top_pa_genes1<-top_pa_genes1%>%
  dplyr::filter(waldStat_2!=0)
topgenes<-rownames(top_pa_genes1[order(top_pa_genes1$waldStat_1,decreasing = T),])[1:50]
pst.ord<-order(m_seu_sce$crv$pseudotime.Lineage1,na.last = NA)
# pseudotime order
heatdata<-heat_count[topgenes,pst.ord]

# draw the heatmap with ComplexHeatmap
# head annotation
anno_col<-fhs_thyrocyte[[c('thca_ident','EMT_score1','TDS_score1')]]
anno_col$ID<-rownames(anno_col)
colnames(anno_col)<-c('thca_ident','EMT_score','TDS_score','ID')
# align the order
# the order of colnames should be the same as the heatdata when cluster_columns=F
anno_col<-anno_col[match(colnames(heatdata),anno_col$ID),]
anno_col<-select(anno_col,c('thca_ident','EMT_score','TDS_score'))
RColorBrewer::brewer.pal(8, "Set2")
col_ha1<-HeatmapAnnotation(df=select(anno_col,c('thca_ident')),
                           EMT_score=anno_lines(select(anno_col,c('EMT_score')),add_points=F,smooth=T,
                                                height = unit(1,'cm'),ylim = c(0,0.06)),
                           TDS_score=anno_lines(select(anno_col,c('TDS_score')),add_points=F,smooth=T,
                                                height = unit(1,'cm'),ylim = c(-0.1,0.1)),
                           show_annotation_name = T,
                           col = list(thca_ident=c('NMB-Thyrocyte'='#f87b71','HLA-Thyrocyte'='#c29c24',
                                                   'NMU-Thyrocyte'='#4db422','IGFBP5-Thyrocyte'='#00bf95',
                                                   'MSMP-Thyrocyte'='#1cb3e8','MTs-Thyrocyte'='#a988fb',
                                                   'RGS5-Thyrocyte'='#fc65d5')),
                           gap = unit(1.25,'mm'))

Heatmap(log1p(as.matrix(heatdata)),
        name = 'expression',
        top_annotation=col_ha1,
        show_column_names = F,
        cluster_columns = F,
        row_km = 4,
        row_title_rot = 0,
        border = T)

# lineage 2
ggplot(slingshot_df, aes(x = slingPseudotime_2, y = thca_ident,
                         colour = thca_ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime 2")

top_pa_genes2<-top_pa_genes
# Replace all NA values with 1 so those genes won't be excluded in the following filtering
top_pa_genes2[is.na(top_pa_genes2)] <- 1
top_pa_genes2<-top_pa_genes2%>%
  dplyr::filter(waldStat_1!=0)
topgenes<-rownames(top_pa_genes2[order(top_pa_genes2$waldStat_2,decreasing = T),])[1:50]
pst.ord<-order(m_seu_sce$crv$pseudotime.Lineage2,na.last = NA)
# pseudotime order
heatdata<-heat_count[topgenes,pst.ord]

# draw the heatmap with ComplexHeatmap
# head annotation

# align the order
# the order of colnames should be the same as the heatdata when cluster_columns=F
anno_col$ID<-rownames(anno_col)
anno_col<-anno_col[match(colnames(heatdata),anno_col$ID),]
anno_col<-select(anno_col,c('thca_ident','EMT_score','TDS_score'))
RColorBrewer::brewer.pal(8, "Set2")
col_ha2<-HeatmapAnnotation(df=select(anno_col,c('thca_ident')),
                           EMT_score=anno_lines(select(anno_col,c('EMT_score')),add_points=F,smooth=T,
                                                height = unit(1,'cm'),ylim = c(0,0.06)),
                           TDS_score=anno_lines(select(anno_col,c('TDS_score')),add_points=F,smooth=T,
                                                height = unit(1,'cm'),ylim = c(-0.1,0.1)),
                           show_annotation_name = T,
                           col = list(thca_ident=c('NMB-Thyrocyte'='#f87b71','HLA-Thyrocyte'='#c29c24',
                                                   'NMU-Thyrocyte'='#4db422','IGFBP5-Thyrocyte'='#00bf95',
                                                   'MSMP-Thyrocyte'='#1cb3e8','MTs-Thyrocyte'='#a988fb',
                                                   'RGS5-Thyrocyte'='#fc65d5')),
                           gap = unit(1.25,'mm'))

# the order of colnames should be the same as the heatdata when cluster_columns=F
Heatmap(log1p(as.matrix(heatdata)),
        name = 'expression',
        top_annotation=col_ha2,
        show_column_names = F,
        cluster_columns = F,
        row_km = 4,
        row_title_rot = 0,
        border = T)

# lineage 3
# in version 2 there is no lineage 3
ggplot(slingshot_df, aes(x = slingPseudotime_3, y = thca_ident,
                         colour = thca_ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime 3")

top_pa_genes3<-top_pa_genes
# Replace all NA values with 1 so those genes won't be excluded in the following filtering
top_pa_genes3[is.na(top_pa_genes3)] <- 1
top_pa_genes3<-top_pa_genes3%>%
  dplyr::filter(waldStat_1!=0&waldStat_2!=0)
topgenes<-rownames(top_pa_genes3[order(top_pa_genes3$waldStat_3,decreasing = T),])[1:50]
pst.ord<-order(m_seu_sce$crv$pseudotime.Lineage3,na.last = NA)
# pseudotime order
heatdata<-heat_count[topgenes,pst.ord]

# draw the heatmap with ComplexHeatmap
# head annotation
# align the order
# the order of colnames should be the same as the heatdata when cluster_columns=F
anno_col$ID<-rownames(anno_col)
anno_col<-anno_col[match(colnames(heatdata),anno_col$ID),]
anno_col<-select(anno_col,c('thca_ident','EMT_score','TDS_score'))
RColorBrewer::brewer.pal(8, "Set2")
col_ha3<-HeatmapAnnotation(df=select(anno_col,c('thca_ident')),
                           EMT_score=anno_lines(select(anno_col,c('EMT_score')),add_points=F,smooth=T,
                                                height = unit(1,'cm'),ylim = c(0,0.06)),
                           TDS_score=anno_lines(select(anno_col,c('TDS_score')),add_points=F,smooth=T,
                                                height = unit(1,'cm'),ylim = c(-0.1,0.2)),
                           show_annotation_name = T,
                           col = list(thca_ident=c('NMB-Thyrocyte'='#66C2A5','HLA-Thyrocyte'='#FC8D62',
                                                   'NMU-Thyrocyte'='#8DA0CB','IGFBP5-Thyrocyte'='#E78AC3',
                                                   'MSMP-Thyrocyte'='#A6D854','MTs-Thyrocyte'='#FFD92F',
                                                   'RGS5-Thyrocyte'='#E5C494')),
                           gap = unit(1.25,'mm'))

# the order of colnames should be the same as the heatdata when cluster_columns=F
Heatmap(log1p(as.matrix(heatdata)),
        name = 'expression',
        top_annotation=col_ha3,
        show_column_names = F,
        cluster_columns = F,
        row_km = 4,
        row_title_rot = 0,
        border = T)
#####

####SCISSOR Analysis####
# the possible link between scRNAseq and bulkseq based on survival data
library(Scissor)
# load the TCGA THCA data, survival data and bulkseq result
# survival_cancer=expr_matrix, clin_meta=survival data
load('TCGA_THCA_phenotype&expr.RData')
# modify the expr_matrix
bulk_dataset<-survival_cancer[,5:810]
bulk_dataset<-t(as.matrix(bulk_dataset))
# modify the phenotype/survival data
phenotype<-clin_meta[,c(3,2)]
colnames(phenotype)<-c('time','status')

# Execute Scissor to select the informative cells
# use Scissor to select the phenotype-associated cell subpopulations, which is fitted by a Cox regression model
# when the five-number summary is <0.01, the result is unreliable
infos1 <- Scissor(bulk_dataset, fhs_thyrocyte, phenotype, alpha = 0.05, 
                  family = "cox", Save_file = 'Scissor_THCA_survival.RData')
# change the Scissor funtion so that it can be applied with SCT normalization
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
infos3 <- Scissor_change(bulk_dataset, fhs_thyrocyte, phenotype, alpha = seq(1,10,2)/1000,cutoff = 0.2, 
                         family = "cox", Load_file = 'Scissor_THCA_survival.RData')
# info3 result: alpha = 0.003, selected cells 26.7%; alpha = 0.005, selected cells 19.7%
# when alpha = 0.05, the %selected cell is 4.88%. Change the alpha to generate a result ~20%
infos1 <- Scissor_change(bulk_dataset, fhs_thyrocyte, phenotype, alpha = 0.005, 
                  family = "cox", Save_file = 'Scissor_THCA_survival.RData')
# Scissor identifies 105 Scissor+ cells associated with worse survival 
# and 996 Scissor- cells associated with good survival

# visualizing the Scissor result in Seurat
Scissor_select <- rep(0, ncol(fhs_thyrocyte))
names(Scissor_select) <- colnames(fhs_thyrocyte)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
fhs_thyrocyte <- AddMetaData(fhs_thyrocyte, metadata = Scissor_select, col.name = "Scissor_OS")
DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_OS', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_OS', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1),
        split.by = 'annotation', ncol = 3)

table(Idents(fhs_thyrocyte),fhs_thyrocyte$Scissor_OS)

rm(fhs_thyrocyte)

## Reliability significance test (heavily time-consuming, 15 h for 50 times)
# setting
# This p-value is less than 0.05, indicating that these associations are reliable. 
numbers <- length(infos1$Scissor_pos) + length(infos1$Scissor_neg)
# load the infos1 result
load('Scissor_THCA_survival.RData')
# in reality nfold should be ~100
result1 <- reliability.test(X, Y, network, alpha = 0.005, family = "cox", cell_num = numbers, n = 10, nfold = 50)
# the result mostly is p=NA & NaN when alpha = 0.05 or 0.005 in infos1

# cell level evaluation
evaluate_summary <- evaluate.cell('Scissor_THCA_survival.RData', infos1, FDR = 0.05, bootstrap_n = 100)
all(evaluate_summary$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary$`Correlation > 0`)) > 50)

# SAVE the result with Scissor result
save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

# the possible link between scRNAseq and bulkseq based on staging data
library(Scissor)
# load the TCGA THCA data, survival data and bulkseq result
# survival_cancer=expr_matrix, clin_meta=survival data
load('TCGA_THCA_phenotype&expr.RData')
# modify the expr_matrix
bulk_dataset<-survival_cancer[,5:810]
bulk_dataset<-t(as.matrix(bulk_dataset))
# modify the phenotype/survival data
phenotype1<-clin_meta[,c(1,7)]
phenotype1$simp_stage <- as.numeric(ifelse(phenotype1$simp_stage %in% c('I', 'II'), 0, 1))
THCA_Stage<-phenotype1$simp_stage
names(THCA_Stage)<-rownames(phenotype1)
tags<-c('early stage','late stage')

# Execute Scissor to select the informative cells
# use Scissor to select the phenotype-associated cell subpopulations, which is fitted by a Cox regression model
# when the five-number summary is <0.01, the result is unreliable
infos2 <- Scissor_change(bulk_dataset, fhs_thyrocyte, THCA_Stage, tag = tags, alpha = 0.2, 
                  family = "binomial", Save_file = 'Scissor_THCA_stage.RData')

Scissor_select <- rep(0, ncol(fhs_thyrocyte))
names(Scissor_select) <- colnames(fhs_thyrocyte)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
fhs_thyrocyte <- AddMetaData(fhs_thyrocyte, metadata = Scissor_select, col.name = "Scissor")
DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1),
        split.by = 'annotation', ncol = 3)

table(Idents(fhs_thyrocyte),fhs_thyrocyte$Scissor)

rm(fhs_thyrocyte)

## Reliability significance test (heavily time-consuming, 15 h for 50 times)
# setting
# This p-value is less than 0.05, indicating that these associations are reliable. 
numbers <- length(infos2$Scissor_pos) + length(infos2$Scissor_neg)
# load the infos1 result
load('Scissor_THCA_stage.RData')
# in reality nfold should be ~100
result2 <- reliability.test(X, Y, network, alpha = 0.2, family = "binomial", cell_num = numbers, n = 10, nfold = 50)
# the result mostly is p=NA & NaN when alpha = 0.2 in infos2

# cell level evaluation
evaluate_summary <- evaluate.cell('Scissor_THCA_stage.RData', infos2, FDR = 0.05, bootstrap_n = 100)
all(evaluate_summary$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary$`Correlation > 0`)) > 50)

# infos3: only IV as late stage
# the possible link between scRNAseq and bulkseq based on staging data
library(Scissor)
# load the TCGA THCA data, survival data and bulkseq result
# survival_cancer=expr_matrix, clin_meta=survival data
load('TCGA_THCA_phenotype&expr.RData')
# modify the expr_matrix
bulk_dataset<-survival_cancer[,5:810]
bulk_dataset<-t(as.matrix(bulk_dataset))
# modify the phenotype/survival data
phenotype2<-clin_meta[,c(1,7)]
phenotype2$simp_stage <- as.numeric(ifelse(phenotype2$simp_stage %in% c('I', 'II','III'), 0, 1))
THCA_Stage<-phenotype2$simp_stage
names(THCA_Stage)<-rownames(phenotype2)
tags<-c('early stage','late stage')

# Execute Scissor to select the informative cells
# use Scissor to select the phenotype-associated cell subpopulations, which is fitted by a Cox regression model
# when the five-number summary is <0.01, the result is unreliable
infos3 <- Scissor_change(bulk_dataset, fhs_thyrocyte, THCA_Stage, tag = tags, alpha = 0.05, 
                         family = "binomial", Save_file = 'Scissor_THCA_stage_v2.RData')

Scissor_select <- rep(0, ncol(fhs_thyrocyte))
names(Scissor_select) <- colnames(fhs_thyrocyte)
Scissor_select[infos3$Scissor_pos] <- 1
Scissor_select[infos3$Scissor_neg] <- 2
fhs_thyrocyte <- AddMetaData(fhs_thyrocyte, metadata = Scissor_select, col.name = "Scissor_stage")
DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_stage', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_stage', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1),
        split.by = 'annotation', ncol = 3)

table(Idents(fhs_thyrocyte),fhs_thyrocyte$Scissor_stage)

rm(fhs_thyrocyte)

## Reliability significance test (heavily time-consuming, 15 h for 50 times)
# setting
# This p-value is less than 0.05, indicating that these associations are reliable. 
numbers <- length(infos3$Scissor_pos) + length(infos3$Scissor_neg)
# load the infos1 result
load('Scissor_THCA_stage_v2.RData')
# in reality nfold should be ~100
result3 <- reliability.test(X, Y, network, alpha = 0.05, family = "binomial", cell_num = numbers, n = 10, nfold = 50)
# 

# cell level evaluation
evaluate_summary_2 <- evaluate.cell('Scissor_THCA_stage_v2.RData', infos3, FDR = 0.05, bootstrap_n = 100)
all(evaluate_summary_2$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary_2$`Correlation > 0`)) > 50)

# SAVE the result with Scissor result
save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

#####abandoned recurrence calculation####
# infos4: recurrence analysis
# the possible link between scRNAseq and bulkseq based on staging data
library(Scissor)
# load the TCGA THCA data, survival data and bulkseq result
# survival_cancer=expr_matrix, clin_meta=survival data
load('TCGA_THCA_phenotype&expr.RData')
# modify the expr_matrix
bulk_dataset<-survival_cancer[,5:810]
bulk_dataset<-t(as.matrix(bulk_dataset))
# modify the phenotype/survival data
phenotype3<-clin_meta[,c(1,12)]
phenotype3$recurrence[is.na(phenotype3$recurrence)] <- 0
phenotype3$recurrence <- as.numeric(ifelse(phenotype3$recurrence %in% c('NO', '0'), 0, 1))
Recurrence<-phenotype3$recurrence
names(Recurrence)<-rownames(phenotype3)
tags<-c('no rec','yes rec')

# Execute Scissor to select the informative cells
# use Scissor to select the phenotype-associated cell subpopulations, which is fitted by a Cox regression model
# when the five-number summary is <0.01, the result is unreliable
Idents(fhs_thyrocyte)<-'annotation'
infos4 <- Scissor_change(bulk_dataset, fhs_thyrocyte, Recurrence, tag = tags, alpha = 0.002, 
                         family = "binomial", Save_file = 'Scissor_THCA_recurrence.RData')

Scissor_select <- rep(0, ncol(fhs_thyrocyte))
names(Scissor_select) <- colnames(fhs_thyrocyte)
Scissor_select[infos4$Scissor_pos] <- 1
Scissor_select[infos4$Scissor_neg] <- 2
fhs_thyrocyte <- AddMetaData(fhs_thyrocyte, metadata = Scissor_select, col.name = "Scissor_recurrence")
DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_recurrence', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))

DimPlot(fhs_thyrocyte, reduction = 'umap', group.by = 'Scissor_recurrence', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1),
        split.by = 'annotation', ncol = 3)

table(Idents(fhs_thyrocyte),fhs_thyrocyte$Scissor_recurrence)

rm(fhs_thyrocyte)

## Reliability significance test (heavily time-consuming, 15 h for 50 times)
# setting
# This p-value is less than 0.05, indicating that these associations are reliable. 
numbers <- length(infos4$Scissor_pos) + length(infos4$Scissor_neg)
# load the infos4 result
load('Scissor_THCA_recurrence.RData')
# in reality nfold should be ~100
result4 <- reliability.test(X, Y, network, alpha = 0.002, family = "binomial", cell_num = numbers, n = 10, nfold = 50)
# did not pass

# cell level evaluation
evaluate_summary_3 <- evaluate.cell('Scissor_THCA_recurrence.RData', infos4, FDR = 0.05, bootstrap_n = 100)
all(evaluate_summary_3$`Mean correlation` & as.numeric(gsub('%','',evaluate_summary_3$`Correlation > 0`)) > 50)

# identify the different-expressed genes between Scissor+ and Scissor- cells
# alpha = 0.05
Idents(fhs_thyrocyte)<-"Scissor_recurrence"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="1",ident.2="2",
                                min.pct = 0.25,logfc.threshold = log(2),min.diff.pct = 0.25)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL
Scissor.de.markers.up<-Scissor.de.markers%>%
  dplyr::filter(Expression=='Up-regulated')

# export the gene list
write.csv(Scissor.de.markers,'Scissor_recurrence_ptc_de_genes.csv')
#####

# SAVE the result with Scissor result
save(fhs_thyrocyte,file="combinedCNV_ptc_TNSC_removed_based_on_SCEVANpred.RData")

####GSEA####
# identify the different-expressed genes between Scissor+ and Scissor- cells
# alpha = 0.05
Idents(fhs_thyrocyte)<-"Scissor_stage"
Scissor.de.markers<-FindMarkers(fhs_thyrocyte, ident.1="1",ident.2="2",
                                min.pct = 0.25,logfc.threshold = log(2),min.diff.pct = 0.25)
Scissor.de.markers$SYMBOL<-rownames(Scissor.de.markers)
# remove those ribosome genes
rb.genes<-rownames(fhs_thyrocyte)[grep("^RP[LS]",rownames(fhs_thyrocyte))]
Scissor.de.markers<-Scissor.de.markers%>%as_tibble()%>%filter(.,!SYMBOL %in% rb.genes)%>%as.data.frame()
rownames(Scissor.de.markers)<-Scissor.de.markers$SYMBOL
Scissor.de.markers.up<-Scissor.de.markers%>%
  dplyr::filter(Expression=='Up-regulated')

# export the gene list
write.csv(Scissor.de.markers,'Scissor_stagev2_ptc_de_genes.csv')
# visualization with volcano plot
Scissor.de.markers<-Scissor.de.markers%>%
  mutate(Expression=case_when(avg_log2FC>=1&p_val_adj<0.05~"Up-regulated",
                              avg_log2FC<=-1&p_val_adj<0.05~"Down-regulated",
                              TRUE~"Unchanged"))
p10<-ggplot(Scissor.de.markers,aes(x=avg_log2FC,y=-log(p_val_adj,10)))+geom_point(aes(color=Expression),size=1)+
  xlab(expression("log"[2]*"FC"))+ylab(expression("-log"[10]*"FDR"))+
  scale_color_manual(values=c("dodgerblue3","gray50","firebrick3"))+
  guides(color=guide_legend(override.aes = list(size=1.5)))
p10
# identify the top 55 genes
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
p11<-p10+geom_label(data=top_Sci_genes,mapping=aes(x=avg_log2FC,y=-log(p_val_adj,10),
                                             label=SYMBOL),size=3)
p11

# GSEA
### enrichment analysis
geneListT<-Scissor.de.markers[,2]
names(geneListT)<-as.character(Scissor.de.markers[,6])
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
emC2_KEGG<-emC2
emC2_KEGG@result<-emC2_KEGG@result[grep("KEGG_", emC2@result$Description),]
enrichplot::dotplot(emC2_KEGG,showCategory = 5,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.1,0.6)

C5BP_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:BP")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5BP<-GSEA(geneListT,TERM2GENE=C5BP_t2g)
enrichplot::dotplot(t.emC5BP,showCategory = 10,split=".sign")+facet_grid(~.sign)

C5CC_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:CC")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5CC<-GSEA(geneListT,TERM2GENE=C5CC_t2g)
enrichplot::dotplot(t.emC5CC,showCategory = 15,split=".sign")+facet_grid(~.sign)
# no term is selected

C5MF_t2g<-msigdbr(species = "Homo sapiens",category="C5",subcategory = "GO:MF")%>%
  dplyr::select(gs_name,gene_symbol)
t.emC5MF<-GSEA(geneListT,TERM2GENE=C5MF_t2g)
enrichplot::dotplot(t.emC5MF,showCategory = 10,split=".sign")+facet_grid(~.sign)
# no term is selected

####DEGs Analysis####
library(DESeq2)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(stringr)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ComplexHeatmap)

# load the data
load(file = 'exprSet for early late stage DEGs.RData')

countData<-exprSet
# gene filtration has been done
# Identify the samples at late stage and early stage
Stage_Ident<-ifelse(clin_meta$simp_stage=='IV','late','early')
Stage_Ident<-as.factor(Stage_Ident)
names(Stage_Ident)<-colnames(countData)
Stage_Ident<-as.data.frame(Stage_Ident)

# remove the samples with NA in simp_stage
Stage_Ident<-na.omit(Stage_Ident)
s = intersect(rownames(Stage_Ident),colnames(countData));length(s)
countData<-countData[,s]

# PCA analysis
t_countData<-as.data.frame(t(countData))
t_countData$tissue=Stage_Ident$Stage_Ident
t_countData_pca<-PCA(t_countData[,-19673],graph=T)
fviz_screeplot(t_countData_pca,addlabels=T)
var<-get_pca_var(t_countData_pca)
fviz_pca_ind(t_countData_pca,
             mean.point=F,#去除分组的中心点
             label="none",
             habillage=as.factor(t_countData$tissue),#根据样本类型着色,should be a factor
             addEllipses = T)#添加边界线

# construct a DESeqDataSet object
# design means the column including treatment (design factor). The design column has to be factor
dds<-DESeqDataSetFromMatrix(countData=countData,colData=Stage_Ident,
                            design=~Stage_Ident)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
# make a logical value keep for the filtering
keep<-rowSums(counts(dds))>=10
dds<-dds[keep,]

# set a factor level. Untreated as reference level
dds$TN_Ident<-relevel(dds$Stage_Ident,ref="early")

# run DESeq
dds<-DESeq(dds)
res<-results(dds,alpha=0.01)

res 

# explore results
summary(res)
# out of 19672 with nonzero total read count
# adjusted p-value < 0.01
# LFC > 0 (up)       : 332, 1.7%
# LFC < 0 (down)     : 808, 4.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 10)

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
# 1 Down-regulated    FDR 0.001   237
# 2 Down-regulated     FDR 0.01   142
# 3 Down-regulated     FDR 0.05    76
# 4      Unchanged    Unchanged 19046
# 5   Up-regulated    FDR 0.001    79
# 6   Up-regulated     FDR 0.01    53
# 7   Up-regulated     FDR 0.05    39

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

p8<-ggplot(res_result,aes(x=log2FoldChange,y=-log(padj,10)))+geom_point(aes(color=Expression),size=2/5)+
  xlab(expression("log"[2]*"FC"))+ylab(expression("-log"[10]*"FDR"))+
  scale_color_manual(values=c("dodgerblue3","gray50","firebrick3"))+
  guides(color=guide_legend(override.aes = list(size=1.5)))
p8

p9<-p8+geom_label(data=top_genes,mapping=aes(x=log2FoldChange,y=-log(padj,10),
                                             label=sample_ID),size=2)
p9

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
emC2_KEGG_mod<-emC2_KEGG
rownames(emC2_KEGG@result)
emC2_KEGG_mod@result<-emC2_KEGG@result[-c(3),]
enrichplot::dotplot(emC2_KEGG_mod,showCategory = 5,split=".sign")+facet_grid(~.sign)+ggplot2::xlim(0.1,0.6)

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
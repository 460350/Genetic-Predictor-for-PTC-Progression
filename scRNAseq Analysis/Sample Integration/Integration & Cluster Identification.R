library(tidyverse)
library(Seurat)
library(glmGamPoi)
library(patchwork)
library(uwot)
library(cowplot)
library(harmony)
library(clustree)

####merge####
f_ptcTN<-merge(x = f_ptc2T,
               y = f_ptc3T,
               project = "thyroidT")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc5T,
               project = "thyroidT")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc10T,
               project = "thyroidT")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc2LeN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc3LeN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc5RiN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc6RiN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc7RiN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc10RiN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc11RiN,
               project = "thyroidN")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc4SC,
               project = "thyroidSC")
f_ptcTN<-merge(x = f_ptcTN,
               y = f_ptc11SC,
               project = "thyroidSC")
####

####normalization and feature selection####
f_ptcTN<-f_ptcTN%>%
  FindVariableFeatures(selection.method = "vst",nfeatures = 2000)%>%
  ScaleData()%>%
  SCTransform(vars.to.regress = c("mitoRatio"), vst.flavor = "v2")

# calculate PCs using variable features determined by SCTransform
f_ptcTN<-RunPCA(f_ptcTN,assay = "SCT",npcs = 50)
####

####integration####
fh_ptcTN<-RunHarmony(f_ptcTN,
                     group.by.vars = "seq_folder",
                     reduction = "pca",assay.use = "SCT",reduction.save = "harmony")
####

####clustering & visualization####
fh_ptcTN<-RunUMAP(object = fh_ptcTN,reduction = "harmony",assay = "SCT",dims = 1:40)%>%
  FindNeighbors(reduction = "harmony")%>%
  FindClusters(resolution = c(0.2,0.4,0.6,0.8))

library(clustree)
clustree::clustree(fh_ptcTN@meta.data,prefix="SCT_snn_res.")

Idents(object = fh_ptcTN)<-"SCT_snn_res.0.2"

# change metadata
fh_ptcTNMetadata<-fh_ptcTN@meta.data

fh_ptcTNMetadata$origin<-fh_ptcTNMetadata$seq_folder
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc2t"]<-"Tumor"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc3t"]<-"Tumor"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc5t"]<-"Tumor"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc10t"]<-"Tumor"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc2len"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc3len"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc5rin"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc6rin"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc7rin"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc10rin"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc11rin"]<-"lymphNode"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc4sc"]<-"subcutaneous"
fh_ptcTNMetadata$origin[fh_ptcTNMetadata$origin == "ptc11sc"]<-"subcutaneous"
#add metatdata back to Seurat object
fh_ptcTN@meta.data<-fh_ptcTNMetadata

fh_ptcTN<-RunUMAP(object = fh_ptcTN,reduction = "harmony",assay = "SCT",dims = 1:40)
DimPlot(fh_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        raster = F)
DimPlot(fh_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        group.by = "origin",
        raster = F)
DimPlot(fh_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        group.by = "seq_folder",
        raster = F)
DimPlot(fh_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        split.by = "seq_folder",
        ncol = 2)
# incorporate the cell cycle score
fh_ptcTN<-CellCycleScoring(fh_ptcTN,
                           g2m.features = g2m_genes,
                           s.features = s_genes)
DimPlot(fh_ptcTN,reduction = "pca",
        label = T,
        split.by = "Phase",
        raster = F)
####

####Overall QC####
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio","riboRatio")
FeaturePlot(fh_ptcTN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)
RidgePlot(fh_ptcTN, features = c("PCNA","TOP2A","MCM6","MKI67"),ncol=2)
## G2M.Score & S.Score enriched in certain subgroup

####removal the cell cycle effect####
fh_ptcTN@meta.data$CC.diff<-fh_ptcTN@meta.data$S.Score-fh_ptcTN@meta.data$G2M.Score
fh_ptcTN<-SCTransform(fh_ptcTN, assay='RNA',new.assay.name='SCT',
                      vars.to.regress=c("CC.diff","mitoRatio","riboRatio"))
ElbowPlot(fh_ptcTN, ndims = 50)
fh_ptcTN<-fh_ptcTN%>%
  RunPCA(assay = 'SCT',npcs = 50)%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.2,0.4,0.6,0.8,1.0))
# check the result
View(fh_ptcTN@meta.data)
clustree::clustree(fh_ptcTN@meta.data,prefix="SCT_snn_res.")
#####

####SAVE####
save(fh_ptcTN,file="integratedData/combined_ptc_TNSC.RData")
####

####annotation####
annotations<-read.csv("annotation.csv")

Idents(object = fh_ptcTN)<-"SCT_snn_res.0.2"

fhptcTN_markers<-FindAllMarkers(object = fh_ptcTN,
                                only.pos = T,
                                logfc.threshold = 0.25)
# Combine marker with gene description
fhptcTN_markers<-fhptcTN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 20 markers per cluster
fhptcTN_top20<-fhptcTN_markers%>%
  group_by(cluster)%>%
  top_n(n = 20, wt = avg_log2FC)

# Plot the specific marker distribution in cluster
FeaturePlot(object = fh_ptcTN, 
            features = c("PTHLH"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE,
            split.by = "origin")
VlnPlot(object = fh_ptcTN,
        features = c("CD3D","CD3E"))

####reset the cluster identity to default####
Idents(object = fh_ptcTN)<-"SCT_snn_res.0.2"

####cluster identity for 0.2####
fh_ptcTN<-RenameIdents(object = fh_ptcTN,
                       "0"="Thyrocyte",
                       "1"="T cell",
                       "2"="Myeloid cell",
                       "3"="T cell",
                       "4"="B cell",
                       "5"="T cell/NK",
                       "6"="Thyrocyte",
                       "7"="Thyrocyte",
                       "8"="Thyrocyte",
                       "9"="Fibroblast",
                       "10"="T cell",
                       "11"="Thyrocyte",
                       "12"="Endothelial cell",
                       "13"="B cell",
                       "14"="T cell/B cell",
                       "15"="Plasma cell",
                       "16"="Fibroblast",
                       "17"="Thyrocyte",
                       "18"="pDC",
                       "19"="Plasma cell",
                       "20"="Thyrocyte")
p1<-DimPlot(fh_ptcTN,
            reduction = "umap",
            label = T,
            label.size = 4,
            repel = T,
            raster = F)
DimPlot(fh_ptcTN,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T,
        split.by = "origin",
        raster = F)
p1

fh_ptcTN@meta.data$overall_ident<-Idents(fh_ptcTN)

ptc_features<-c("TG","TSHR","TPO","TRIP6","KRT19","CD3D","CD3E","LYZ","CST3","CD79A","MS4A1",
                "IGKC","IGHG3","MZB1","NKG7","GNLY","COL1A2","COL6A2","VWF","ESAM","PECAM1","FCER1G")

p2<-DotPlot(fh_ptcTN,features = ptc_features,cols = c("blue","red"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()
p2

VlnPlot(fh_ptcTN,features = c("TG"),raster = F)

DoHeatmap(subset(fh_ptcTN,downsample=100),features = ptc_features,
          cells = 1:500, size = 4, angle = 90)+NoLegend()+
  scale_fill_gradientn(colors = c("blue", "white", "red"))

## bar plot
# calculation
cellRatio<-prop.table(table(Idents(fh_ptcTN),fh_ptcTN$orig.ident),margin = 2)%>%
  as.data.frame()
cellRatio
colorCount<-length(unique(cellRatio$Var1))
ggplot(cellRatio)+
  geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat = "identity",width = 0.7,size=0.5,colour='#222222')+
  theme_classic()+
  labs(x='Sample',y='Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color = "black",size=0.5,linetype = "solid"))
####

####SAVE####
save(fh_ptcTN,file="integratedData/combined_ptc_TNSC.RData")
####
library(Seurat)
library(tidyverse)
library(clustree)
library(SCEVAN)

####before split####
# Check the percentage
summary(as.factor(fh_ptcTN@meta.data$orig.ident))
prop.table(table(Idents(fh_ptcTN)))
####

####split####
# Immune cell
fh_Immune<-subset(x=fh_ptcTN,
                  idents=c("Thyrocyte","Fibroblast","Endothelial cell"),
                  invert=T)

DimPlot(fh_Immune,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)

FeaturePlot(object = fh_Immune, 
            features = c("KRT18"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

save(fh_Immune,file="integratedData/combined_ptcImmune_cell.RData")

# nonImmune cell
fh_nonImmune<-subset(x=fh_ptcTN,
                  idents=c("Thyrocyte","Fibroblast","Endothelial cell"))

DimPlot(fh_nonImmune,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)

FeaturePlot(object = fh_nonImmune, 
            features = c("CD3E","CD3D"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

FeaturePlot(object = fh_nonImmune, 
            features = c("TG","KRT18"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

####SAVE####
save(fh_nonImmune,file="integratedData/combined_ptcnonImmune_cell.RData")
####

####further trimming####
# certain group of cells need to be removed 
fh_nonImmune<-SCTransform(fh_nonImmune, assay='RNA',new.assay.name='SCT',
                      vars.to.regress=c("CC.diff","mitoRatio","riboRatio"))
ElbowPlot(fh_nonImmune, ndims = 50)
fh_nonImmune<-fh_nonImmune%>%
  RunPCA(assay = 'SCT',npcs = 50)%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.2,0.4,0.6))
# check the result
View(fh_nonImmune@meta.data)
clustree::clustree(fh_nonImmune@meta.data,prefix="SCT_snn_res.")

Idents(object = fh_nonImmune)<-"SCT_snn_res.0.2"
DimPlot(fh_nonImmune,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)

# isolate cluster 4 for immune cell
immuneCell<-subset(x=fh_nonImmune,idents=c("4"))
save(immuneCell,file="integratedData/combined_ptcNonImmune_ICsub.RData")

# real nonImmune clusters
fh_NonImmune<-subset(x=fh_nonImmune,idents=c("4"),
                     invert=T)
Idents(fh_NonImmune)<-fh_NonImmune@meta.data$overall_ident
DimPlot(fh_NonImmune,
        reduction = "umap",
        label = T,
        label.size = 4,
        repel = T)
FeaturePlot(object = fh_NonImmune, 
            features = c("CD3E","CD3D"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

save(fh_NonImmune,file="integratedData/combined_ptcnonImmune_cell.RData")
####

####SPLIT####
summary(as.factor(fh_NonImmune@meta.data$orig.ident))

###ptc10rin   ptc10t ptc11rin  ptc11sc    ptc2t  ptc3len    ptc3t   ptc4sc  ptc5rin    ptc5t 
### 5642     9989     1552     5318        1       70      511      245     6010     3850 
###ptc6rin  ptc7rin 
###6326     3157 

# split them based on original idents:ptc2t+ptc3len, ptc3t, ptc5rin, ptc5t, ptc6rin
Idents(object = fh_NonImmune)<-"orig.ident"
ptc2_3tn<-subset(fh_NonImmune,idents=c('ptc2t','ptc3len'))
ptc3t<-subset(fh_NonImmune,idents= 'ptc3t')
ptc4sc<-subset(fh_NonImmune,idents= 'ptc4sc')
ptc5rin<-subset(fh_NonImmune,idents= 'ptc5rin')
ptc5t<-subset(fh_NonImmune,idents= 'ptc5t')
ptc6rin<-subset(fh_NonImmune,idents= 'ptc6rin')
ptc7rin<-subset(fh_NonImmune,idents= 'ptc7rin')
ptc10rin<-subset(fh_NonImmune,idents= 'ptc10rin')
ptc10t<-subset(fh_NonImmune,idents= 'ptc10t')

####ptc2t+3len###
ptc2_3tn_mtx<-ptc2_3tn@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc2_3tn.results <- SCEVAN::pipelineCNA(ptc2_3tn_mtx, sample = "ptc2t_3len",
                                       par_cores = 10,
                                       SUBCLONES = F, plotTree = F)

View(ptc2_3tn.results)

ptc2_3tn<-AddMetaData(ptc2_3tn,metadata = ptc2_3tn.results,
                     col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc2_3tn,file="Split Sample/ptc2t3len.RData")
rm(ptc2_3tn)
####

####ptc3t####
ptc3t_mtx<-ptc3t@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc3t.results <- SCEVAN::pipelineCNA(ptc3t_mtx, sample = "ptc3t",
                                     par_cores = 10,
                                     SUBCLONES = F, plotTree = F)

View(ptc3t.results)

ptc3t<-AddMetaData(ptc3t,metadata = ptc3t.results[,-3],
                   col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc3t,file="Split Sample/ptc3t.RData")
rm(ptc3t)
####
####ptc4sc####
ptc4sc_mtx<-ptc4sc@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc4sc.results <- SCEVAN::pipelineCNA(ptc4sc_mtx, sample = "ptc4sc",
                                     par_cores = 10,
                                     SUBCLONES = F, plotTree = F)

View(ptc4sc.results)

ptc4sc<-AddMetaData(ptc4sc,metadata = ptc4sc.results,
                   col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc4sc,file="Split Sample/ptc4sc.RData")
rm(ptc4sc)
####
####ptc5t####
ptc5t_mtx<-ptc5t@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc5t.results <- SCEVAN::pipelineCNA(ptc5t_mtx, sample = "ptc5t",
                                     par_cores = 10,
                                     SUBCLONES = F, plotTree = F)

View(ptc5t.results)

ptc5t<-AddMetaData(ptc5t,metadata = ptc5t.results,
                   col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc5t,file="Split Sample/ptc5t.RData")
rm(ptc5t)
####
####ptc5rin####
ptc5rin_mtx<-ptc5rin@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc5rin.results <- SCEVAN::pipelineCNA(ptc5rin_mtx, sample = "ptc5rin",
                                     par_cores = 10,
                                     SUBCLONES = F, plotTree = F)

View(ptc5rin.results)

ptc5rin<-AddMetaData(ptc5rin,metadata = ptc5rin.results,
                   col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc5rin,file="Split Sample/ptc5rin.RData")
rm(ptc5rin)
####
####ptc6rin####
ptc6rin_mtx<-ptc6rin@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc6rin.results <- SCEVAN::pipelineCNA(ptc6rin_mtx, sample = "ptc6rin",
                                       par_cores = 10,
                                       SUBCLONES = F, plotTree = F)

View(ptc6rin.results)

ptc6rin<-AddMetaData(ptc6rin,metadata = ptc6rin.results,
                     col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc6rin,file="Split Sample/ptc6rin.RData")
rm(ptc6rin)
####
####ptc7rin####
ptc7rin_mtx<-ptc7rin@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc7rin.results <- SCEVAN::pipelineCNA(ptc7rin_mtx, sample = "ptc7rin",
                                       par_cores = 10,
                                       SUBCLONES = F, plotTree = F)

View(ptc7rin.results)

ptc7rin<-AddMetaData(ptc7rin,metadata = ptc7rin.results,
                     col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc7rin,file="Split Sample/ptc7rin.RData")
rm(ptc7rin)
####
####ptc10rin####
ptc10rin_mtx<-ptc10rin@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc10rin.results <- SCEVAN::pipelineCNA(ptc10rin_mtx, sample = "ptc10rin",
                                       par_cores = 10,
                                       SUBCLONES = F, plotTree = F)

View(ptc10rin.results)

ptc10rin<-AddMetaData(ptc10rin,metadata = ptc10rin.results,
                     col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc10rin,file="Split Sample/ptc10rin.RData")
rm(ptc10rin)
####
####ptc10t####
ptc10t_mtx<-ptc10t@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc10t.results <- SCEVAN::pipelineCNA(ptc10t_mtx, sample = "ptc10t",
                                     par_cores = 10,
                                     SUBCLONES = F, plotTree = F)

View(ptc10t.results)

ptc10t<-AddMetaData(ptc10t,metadata = ptc10t.results,
                   col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc10t,file="Split Sample/ptc10t.RData")
rm(ptc10t)
####
####
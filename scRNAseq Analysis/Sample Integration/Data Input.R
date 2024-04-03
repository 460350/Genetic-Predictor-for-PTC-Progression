####update note####
# 2023-5-5, add the code to calculate the riboratio and update the processed files.
#      only update the tumor, lymph node, and subcutaneous samples.
#####

library(tidyverse)
library(Seurat)
library(patchwork)
library(uwot)
library(cowplot)

#raw sample data are from 10x cellranger

###########Thyroid 1 Tumor################
#load tumor-pre dataset
ptc1T.data<-Read10X(data.dir = "scRNAseq Data/PTC1+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc1T<-CreateSeuratObject(counts = ptc1T.data,project = "ptc1t", 
                               min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc1T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc1T$log10GenesPerUMI <- log10(ptc1T$nFeature_RNA) / log10(ptc1T$nCount_RNA)
# Compute percent mito ratio
ptc1T$mitoRatio <- PercentageFeatureSet(object = ptc1T, pattern = "^MT-")
ptc1T$mitoRatio <- ptc1T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc1T$riboRatio <- PercentageFeatureSet(object = ptc1T, pattern = "^RP[LS]")
ptc1T$riboRatio <- ptc1T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc1TMetadata<-ptc1T@meta.data
#add cell IDs to metadata
ptc1TMetadata$cells<-rownames(ptc1TMetadata)
#rename columns
ptc1TMetadata<-ptc1TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc1T@meta.data<-ptc1TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc1TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc1TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc1TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc1TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc1TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc1T <- subset(x = ptc1T, 
                       subset= (nUMI >= 500) & 
                         (nGene >= 250) & 
                         (log10GenesPerUMI > 0.80) & 
                         (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC1tCounts<-GetAssayData(object = f_ptc1T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC1tNonzero <- tPTC1tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc1keep_genes <- Matrix::rowSums(tPTC1tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC1tCounts <- tPTC1tCounts[ptc1keep_genes, ]
# Reassign to filtered Seurat object
f_ptc1T <- CreateSeuratObject(f_tPTC1tCounts, meta.data = f_ptc1T@meta.data,
                              project = "ptc1t")
# RE_ASSESS QC METRICS
f_ptc1T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc1T,file="processedData/filtered_PTC1_T.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc1T_phase<-NormalizeData(f_ptc1T)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc1T_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc1T<-SCTransform(f_ptc1T,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc1T<-f_ptc1T%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc1T@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc1T)<-"SCT_snn_res.0.6"
f_ptc1T<-f_ptc1T%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc1T,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc1T<-CellCycleScoring(f_ptc1T,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc1T,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc1T,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:15),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc1T_PCdata<-FetchData(f_ptc1T,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp1t_umap_label <- FetchData(f_ptc1T,
                            vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:15), function(pc){
  ggplot(f_ptc1T_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp1t_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc1T[["pca"]], dims = 1:15, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp1t_markers<-FindAllMarkers(object = f_ptc1T,
                            only.pos = T,
                            logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp1t_markers<-fp1t_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp1t_top10<-fp1t_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
## Unsupervised Cluster Identification


###########Tumor 1 Normal Tis/Para Tum###############
#load para tumor/normal tissue dataset
ptc1P.data<-Read10X(data.dir = "scRNAseq Data/PTC1+P/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc1P<-CreateSeuratObject(counts = ptc1P.data,project = "ptc1p", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc1P@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc1P$log10GenesPerUMI <- log10(ptc1P$nFeature_RNA) / log10(ptc1P$nCount_RNA)
# Compute percent mito ratio
ptc1P$mitoRatio <- PercentageFeatureSet(object = ptc1P, pattern = "^MT-")
ptc1P$mitoRatio <- ptc1P@meta.data$mitoRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc1PMetadata<-ptc1P@meta.data
#add cell IDs to metadata
ptc1PMetadata$cells<-rownames(ptc1PMetadata)
#rename columns
ptc1PMetadata<-ptc1PMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc1P@meta.data<-ptc1PMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc1PMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc1PMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc1PMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc1PMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc1PMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc1P <- subset(x = ptc1P, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC1pCounts<-GetAssayData(object = f_ptc1P,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC1pNonzero <- tPTC1pCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc1pkeep_genes <- Matrix::rowSums(tPTC1pNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC1pCounts <- tPTC1pCounts[ptc1pkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc1P <- CreateSeuratObject(f_tPTC1pCounts, meta.data = f_ptc1P@meta.data,
                              project = "ptc1p")
# RE_ASSESS QC METRICS
f_ptc1P@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc1P,file="processedData/filtered_PTC1_P.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc1P_phase<-NormalizeData(f_ptc1P)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc1P_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc1P<-SCTransform(f_ptc1P,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc1P<-f_ptc1P%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc1P@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc1P)<-"SCT_snn_res.0.6"
f_ptc1P<-f_ptc1P%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc1P,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc1P<-CellCycleScoring(f_ptc1P,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc1P,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc1P,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc1P_PCdata<-FetchData(f_ptc1P,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp1p_umap_label <- FetchData(f_ptc1P,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc1P_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp1p_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc1P[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp1p_markers<-FindAllMarkers(object = f_ptc1P,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp1p_markers<-fp1p_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp1p_top10<-fp1p_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification


###########Thyroid 2 Tumor#################
#load tumor-pre dataset
ptc2T.data<-Read10X(data.dir = "scRNAseq Data/PTC2+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc2T<-CreateSeuratObject(counts = ptc2T.data,project = "ptc2t", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc2T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc2T$log10GenesPerUMI <- log10(ptc2T$nFeature_RNA) / log10(ptc2T$nCount_RNA)
# Compute percent mito ratio
ptc2T$mitoRatio <- PercentageFeatureSet(object = ptc2T, pattern = "^MT-")
ptc2T$mitoRatio <- ptc2T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc2T$riboRatio <- PercentageFeatureSet(object = ptc2T, pattern = "^RP[LS]")
ptc2T$riboRatio <- ptc2T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc2TMetadata<-ptc2T@meta.data
#add cell IDs to metadata
ptc2TMetadata$cells<-rownames(ptc2TMetadata)
#rename columns
ptc2TMetadata<-ptc2TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc2T@meta.data<-ptc2TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc2TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc2TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc2TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc2TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc2TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc2T <- subset(x = ptc2T, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC2tCounts<-GetAssayData(object = f_ptc2T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC2tNonzero <- tPTC2tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc2tkeep_genes <- Matrix::rowSums(tPTC2tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC2tCounts <- tPTC2tCounts[ptc2tkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc2T <- CreateSeuratObject(f_tPTC2tCounts, meta.data = f_ptc2T@meta.data,
                              project = "ptc2t")
# RE_ASSESS QC METRICS
f_ptc2T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc2T,file="processedData/filtered_PTC2_T.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc2T_phase<-NormalizeData(f_ptc2T)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc2T_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc2T<-SCTransform(f_ptc2T,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc2T<-f_ptc2T%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc2T@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc2T)<-"SCT_snn_res.0.6"
f_ptc2T<-f_ptc2T%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc2T,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc2T<-CellCycleScoring(f_ptc2T,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc2T,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc2T,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:13),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc2T_PCdata<-FetchData(f_ptc2T,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp2t_umap_label <- FetchData(f_ptc2T,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:13), function(pc){
  ggplot(f_ptc2T_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp2t_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc2T[["pca"]], dims = 1:13, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp2t_markers<-FindAllMarkers(object = f_ptc2T,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp2t_markers<-fp2t_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp2t_top10<-fp2t_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
## Unsupervised Cluster Identification


###########Thyroid 2 Normal Tis/Para Tum#############
#load para tumor/normal tissue dataset
ptc2P.data<-Read10X(data.dir = "scRNAseq Data/PTC2+P/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc2P<-CreateSeuratObject(counts = ptc2P.data,project = "ptc2p", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc2P@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc2P$log10GenesPerUMI <- log10(ptc2P$nFeature_RNA) / log10(ptc2P$nCount_RNA)
# Compute percent mito ratio
ptc2P$mitoRatio <- PercentageFeatureSet(object = ptc2P, pattern = "^MT-")
ptc2P$mitoRatio <- ptc2P@meta.data$mitoRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc2PMetadata<-ptc2P@meta.data
#add cell IDs to metadata
ptc2PMetadata$cells<-rownames(ptc2PMetadata)
#rename columns
ptc2PMetadata<-ptc2PMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc2P@meta.data<-ptc2PMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc2PMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc2PMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc2PMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc2PMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc2PMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc2P <- subset(x = ptc2P, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC2pCounts<-GetAssayData(object = f_ptc2P,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC2pNonzero <- tPTC2pCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc2pkeep_genes <- Matrix::rowSums(tPTC2pNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC2pCounts <- tPTC2pCounts[ptc2pkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc2P <- CreateSeuratObject(f_tPTC2pCounts, meta.data = f_ptc2P@meta.data,
                              project = "ptc2p")
# RE_ASSESS QC METRICS
f_ptc2P@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc2P,file="processedData/filtered_PTC2_P.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc2P_phase<-NormalizeData(f_ptc2P)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc2P_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc2P<-SCTransform(f_ptc2P,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc2P<-f_ptc2P%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc2P@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc2P)<-"SCT_snn_res.0.6"
f_ptc2P<-f_ptc2P%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc2P,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc2P<-CellCycleScoring(f_ptc2P,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc2P,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc2P,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:15),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc2P_PCdata<-FetchData(f_ptc2P,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp2p_umap_label <- FetchData(f_ptc2P,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:15), function(pc){
  ggplot(f_ptc2P_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp2p_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc2P[["pca"]], dims = 1:15, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp2p_markers<-FindAllMarkers(object = f_ptc2P,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp2p_markers<-fp2p_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp2p_top10<-fp2p_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification


###########Thyroid 3 Tumor################
#load tumor-pre dataset
ptc3T.data<-Read10X(data.dir = "scRNAseq Data/PTC3+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc3T<-CreateSeuratObject(counts = ptc3T.data,project = "ptc3t", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc3T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc3T$log10GenesPerUMI <- log10(ptc3T$nFeature_RNA) / log10(ptc3T$nCount_RNA)
# Compute percent mito ratio
ptc3T$mitoRatio <- PercentageFeatureSet(object = ptc3T, pattern = "^MT-")
ptc3T$mitoRatio <- ptc3T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc3T$riboRatio <- PercentageFeatureSet(object = ptc3T, pattern = "^RP[LS]")
ptc3T$riboRatio <- ptc3T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc3TMetadata<-ptc3T@meta.data
#add cell IDs to metadata
ptc3TMetadata$cells<-rownames(ptc3TMetadata)
#rename columns
ptc3TMetadata<-ptc3TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc3T@meta.data<-ptc3TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc3TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc3TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc3TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc3TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc3TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc3T <- subset(x = ptc3T, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC3tCounts<-GetAssayData(object = f_ptc3T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC3tNonzero <- tPTC3tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc3tkeep_genes <- Matrix::rowSums(tPTC3tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC3tCounts <- tPTC3tCounts[ptc3tkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc3T <- CreateSeuratObject(f_tPTC3tCounts, meta.data = f_ptc3T@meta.data,
                              project = "ptc3t")
# RE_ASSESS QC METRICS
f_ptc3T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc3T,file="processedData/filtered_PTC3_T.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc3T_phase<-NormalizeData(f_ptc3T)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc3T_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc3T<-SCTransform(f_ptc3T,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc3T<-f_ptc3T%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc3T@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc3T)<-"SCT_snn_res.0.6"
f_ptc3T<-f_ptc3T%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc3T,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc3T<-CellCycleScoring(f_ptc3T,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc3T,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc3T,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:21),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc3T_PCdata<-FetchData(f_ptc3T,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp3t_umap_label <- FetchData(f_ptc3T,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:21), function(pc){
  ggplot(f_ptc3T_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp3t_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc3T[["pca"]], dims = 1:21, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp3t_markers<-FindAllMarkers(object = f_ptc3T,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp3t_markers<-fp3t_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp3t_top10<-fp3t_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
## Unsupervised Cluster Identification


###########Thyroid 3 Normal Tis/Para Tum#########
#load para tumor/normal tissue dataset
ptc3P.data<-Read10X(data.dir = "scRNAseq Data/PTC3+P/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc3P<-CreateSeuratObject(counts = ptc3P.data,project = "ptc3p", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc3P@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc3P$log10GenesPerUMI <- log10(ptc3P$nFeature_RNA) / log10(ptc3P$nCount_RNA)
# Compute percent mito ratio
ptc3P$mitoRatio <- PercentageFeatureSet(object = ptc3P, pattern = "^MT-")
ptc3P$mitoRatio <- ptc3P@meta.data$mitoRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc3PMetadata<-ptc3P@meta.data
#add cell IDs to metadata
ptc3PMetadata$cells<-rownames(ptc3PMetadata)
#rename columns
ptc3PMetadata<-ptc3PMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc3P@meta.data<-ptc3PMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc3PMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc3PMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc3PMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc3PMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc3PMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc3P <- subset(x = ptc3P, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC3pCounts<-GetAssayData(object = f_ptc3P,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC3pNonzero <- tPTC3pCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc3pkeep_genes <- Matrix::rowSums(tPTC3pNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC3pCounts <- tPTC3pCounts[ptc3pkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc3P <- CreateSeuratObject(f_tPTC3pCounts, meta.data = f_ptc3P@meta.data,
                              project = "ptc3p")
# RE_ASSESS QC METRICS
f_ptc3P@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc3P,file="processedData/filtered_PTC3_P.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc3P_phase<-NormalizeData(f_ptc3P)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc3P_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc3P<-SCTransform(f_ptc3P,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc3P<-f_ptc3P%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc3P@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc3P)<-"SCT_snn_res.0.6"
f_ptc3P<-f_ptc3P%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc3P,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc3P<-CellCycleScoring(f_ptc3P,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc3P,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc3P,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:19),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc3P_PCdata<-FetchData(f_ptc3P,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp3p_umap_label <- FetchData(f_ptc3P,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:19), function(pc){
  ggplot(f_ptc3P_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp3p_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc3P[["pca"]], dims = 1:19, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp3p_markers<-FindAllMarkers(object = f_ptc3P,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp3p_markers<-fp3p_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp3p_top10<-fp3p_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification


###########Thyroid 5 Tumor###############
#load tumor-pre dataset
ptc5T.data<-Read10X(data.dir = "scRNAseq Data/PTC5+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc5T<-CreateSeuratObject(counts = ptc5T.data,project = "ptc5t", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc5T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc5T$log10GenesPerUMI <- log10(ptc5T$nFeature_RNA) / log10(ptc5T$nCount_RNA)
# Compute percent mito ratio
ptc5T$mitoRatio <- PercentageFeatureSet(object = ptc5T, pattern = "^MT-")
ptc5T$mitoRatio <- ptc5T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc5T$riboRatio <- PercentageFeatureSet(object = ptc5T, pattern = "^RP[LS]")
ptc5T$riboRatio <- ptc5T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc5TMetadata<-ptc5T@meta.data
#add cell IDs to metadata
ptc5TMetadata$cells<-rownames(ptc5TMetadata)
#rename columns
ptc5TMetadata<-ptc5TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc5T@meta.data<-ptc5TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc5TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc5TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc5TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc5TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc5TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc5T <- subset(x = ptc5T, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC5tCounts<-GetAssayData(object = f_ptc5T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC5tNonzero <- tPTC5tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc5tkeep_genes <- Matrix::rowSums(tPTC5tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC5tCounts <- tPTC5tCounts[ptc5tkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc5T <- CreateSeuratObject(f_tPTC5tCounts, meta.data = f_ptc5T@meta.data,
                              project = "ptc5t")
# RE_ASSESS QC METRICS
f_ptc5T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc5T,file="processedData/filtered_PTC5_T.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc5T_phase<-NormalizeData(f_ptc5T)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc5T_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc5T<-SCTransform(f_ptc5T,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc5T<-f_ptc5T%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc5T@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc5T)<-"SCT_snn_res.0.6"
f_ptc5T<-f_ptc5T%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc5T,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc5T<-CellCycleScoring(f_ptc5T,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc5T,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc5T,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:13),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc5T_PCdata<-FetchData(f_ptc5T,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp5t_umap_label <- FetchData(f_ptc5T,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:13), function(pc){
  ggplot(f_ptc5T_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp5t_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc5T[["pca"]], dims = 1:13, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp5t_markers<-FindAllMarkers(object = f_ptc5T,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp5t_markers<-fp5t_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp5t_top10<-fp5t_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
## Unsupervised Cluster Identification

###########Thyroid 5 Normal Tis/Para Tum###########
#load para tumor/normal tissue dataset
ptc5P.data<-Read10X(data.dir = "scRNAseq Data/PTC5+P/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc5P<-CreateSeuratObject(counts = ptc5P.data,project = "ptc5p", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc5P@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc5P$log10GenesPerUMI <- log10(ptc5P$nFeature_RNA) / log10(ptc5P$nCount_RNA)
# Compute percent mito ratio
ptc5P$mitoRatio <- PercentageFeatureSet(object = ptc5P, pattern = "^MT-")
ptc5P$mitoRatio <- ptc5P@meta.data$mitoRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc5PMetadata<-ptc5P@meta.data
#add cell IDs to metadata
ptc5PMetadata$cells<-rownames(ptc5PMetadata)
#rename columns
ptc5PMetadata<-ptc5PMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc5P@meta.data<-ptc5PMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc5PMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc5PMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc5PMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc5PMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc5PMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc5P <- subset(x = ptc5P, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC5pCounts<-GetAssayData(object = f_ptc5P,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC5pNonzero <- tPTC5pCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc5pkeep_genes <- Matrix::rowSums(tPTC5pNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC5pCounts <- tPTC5pCounts[ptc5pkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc5P <- CreateSeuratObject(f_tPTC5pCounts, meta.data = f_ptc5P@meta.data,
                              project = "ptc5p")
# RE_ASSESS QC METRICS
f_ptc5P@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc5P,file="processedData/filtered_PTC5_P.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc5P_phase<-NormalizeData(f_ptc5P)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc5P_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc5P<-SCTransform(f_ptc5P,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc5P<-f_ptc5P%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc5P@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc5P)<-"SCT_snn_res.0.6"
f_ptc5P<-f_ptc5P%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc5P,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc5P<-CellCycleScoring(f_ptc5P,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc5P,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc5P,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:16),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc5P_PCdata<-FetchData(f_ptc5P,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp5p_umap_label <- FetchData(f_ptc5P,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:16), function(pc){
  ggplot(f_ptc5P_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp5p_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc5P[["pca"]], dims = 1:16, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp5p_markers<-FindAllMarkers(object = f_ptc5P,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp5p_markers<-fp5p_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp5p_top10<-fp5p_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification

###########Thyroid 2 Left N##################
#load LN mets tissue dataset
ptc2LeN.data<-Read10X(data.dir = "scRNAseq Data/PTC2+LeftN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc2LeN<-CreateSeuratObject(counts = ptc2LeN.data,project = "ptc2len", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc2LeN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc2LeN$log10GenesPerUMI <- log10(ptc2LeN$nFeature_RNA) / log10(ptc2LeN$nCount_RNA)
# Compute percent mito ratio
ptc2LeN$mitoRatio <- PercentageFeatureSet(object = ptc2LeN, pattern = "^MT-")
ptc2LeN$mitoRatio <- ptc2LeN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc2LeN$riboRatio <- PercentageFeatureSet(object = ptc2LeN, pattern = "^RP[LS]")
ptc2LeN$riboRatio <- ptc2LeN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc2LeNMetadata<-ptc2LeN@meta.data
#add cell IDs to metadata
ptc2LeNMetadata$cells<-rownames(ptc2LeNMetadata)
#rename columns
ptc2LeNMetadata<-ptc2LeNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc2LeN@meta.data<-ptc2LeNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc2LeNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc2LeNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc2LeNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc2LeNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc2LeNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc2LeN <- subset(x = ptc2LeN, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC2LeNCounts<-GetAssayData(object = f_ptc2LeN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC2LeNNonzero <- tPTC2LeNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc2LeNkeep_genes <- Matrix::rowSums(tPTC2LeNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC2LeNCounts <- tPTC2LeNCounts[ptc2LeNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc2LeN <- CreateSeuratObject(f_tPTC2LeNCounts, meta.data = f_ptc2LeN@meta.data,
                              project = "ptc2len")
# RE_ASSESS QC METRICS
f_ptc2LeN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc2LeN,file="processedData/filtered_PTC2_LeftN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc2LeN_phase<-NormalizeData(f_ptc2LeN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc2LeN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc2LeN<-SCTransform(f_ptc2LeN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc2LeN<-f_ptc2LeN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc2LeN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc2LeN)<-"SCT_snn_res.0.6"
f_ptc2LeN<-f_ptc2LeN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc2LeN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc2LeN<-CellCycleScoring(f_ptc2LeN,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc2LeN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc2LeN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:16),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc2LeN_PCdata<-FetchData(f_ptc2LeN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc2LeN_umap_label <- FetchData(f_ptc2LeN,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:16), function(pc){
  ggplot(f_ptc2LeN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc2LeN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc2LeN[["pca"]], dims = 1:16, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc2LeN_markers<-FindAllMarkers(object = f_ptc2LeN,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc2LeN_markers<-fptc2LeN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc2LeN_top10<-fptc2LeN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification


###########Thyroid 3 Left N###################
#load LN mets tissue dataset
ptc3LeN.data<-Read10X(data.dir = "scRNAseq Data/PTC3+LeftN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc3LeN<-CreateSeuratObject(counts = ptc3LeN.data,project = "ptc3len", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc3LeN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc3LeN$log10GenesPerUMI <- log10(ptc3LeN$nFeature_RNA) / log10(ptc3LeN$nCount_RNA)
# Compute percent mito ratio
ptc3LeN$mitoRatio <- PercentageFeatureSet(object = ptc3LeN, pattern = "^MT-")
ptc3LeN$mitoRatio <- ptc3LeN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc3LeN$riboRatio <- PercentageFeatureSet(object = ptc3LeN, pattern = "^RP[LS]")
ptc3LeN$riboRatio <- ptc3LeN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc3LeNMetadata<-ptc3LeN@meta.data
#add cell IDs to metadata
ptc3LeNMetadata$cells<-rownames(ptc3LeNMetadata)
#rename columns
ptc3LeNMetadata<-ptc3LeNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc3LeN@meta.data<-ptc3LeNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc3LeNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc3LeNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc3LeNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc3LeNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc3LeNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc3LeN <- subset(x = ptc3LeN, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC3LeNCounts<-GetAssayData(object = f_ptc3LeN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC3LeNNonzero <- tPTC3LeNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc3LeNkeep_genes <- Matrix::rowSums(tPTC3LeNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC3LeNCounts <- tPTC3LeNCounts[ptc3LeNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc3LeN <- CreateSeuratObject(f_tPTC3LeNCounts, meta.data = f_ptc3LeN@meta.data,
                                project = "ptc3len")
# RE_ASSESS QC METRICS
f_ptc3LeN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc3LeN,file="processedData/filtered_PTC3_LeftN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc3LeN_phase<-NormalizeData(f_ptc3LeN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc3LeN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc3LeN<-SCTransform(f_ptc3LeN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc3LeN<-f_ptc3LeN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc3LeN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc3LeN)<-"SCT_snn_res.0.6"
f_ptc3LeN<-f_ptc3LeN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc3LeN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc3LeN<-CellCycleScoring(f_ptc3LeN,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc3LeN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc3LeN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:16),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc3LeN_PCdata<-FetchData(f_ptc3LeN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc3LeN_umap_label <- FetchData(f_ptc3LeN,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:16), function(pc){
  ggplot(f_ptc3LeN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc3LeN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc3LeN[["pca"]], dims = 1:16, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc3LeN_markers<-FindAllMarkers(object = f_ptc3LeN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc3LeN_markers<-fptc3LeN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc3LeN_top10<-fptc3LeN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification


###########Thyroid 5 Right N##################
#load LN mets tissue dataset
ptc5RiN.data<-Read10X(data.dir = "scRNAseq Data/PTC5+RightN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc5RiN<-CreateSeuratObject(counts = ptc5RiN.data,project = "ptc5rin", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc5RiN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc5RiN$log10GenesPerUMI <- log10(ptc5RiN$nFeature_RNA) / log10(ptc5RiN$nCount_RNA)
# Compute percent mito ratio
ptc5RiN$mitoRatio <- PercentageFeatureSet(object = ptc5RiN, pattern = "^MT-")
ptc5RiN$mitoRatio <- ptc5RiN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc5RiN$riboRatio <- PercentageFeatureSet(object = ptc5RiN, pattern = "^RP[LS]")
ptc5RiN$riboRatio <- ptc5RiN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc5RiNMetadata<-ptc5RiN@meta.data
#add cell IDs to metadata
ptc5RiNMetadata$cells<-rownames(ptc5RiNMetadata)
#rename columns
ptc5RiNMetadata<-ptc5RiNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc5RiN@meta.data<-ptc5RiNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc5RiNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc5RiNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc5RiNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc5RiNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc5RiNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc5RiN <- subset(x = ptc5RiN, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC5RiNCounts<-GetAssayData(object = f_ptc5RiN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC5RiNNonzero <- tPTC5RiNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc5RiNkeep_genes <- Matrix::rowSums(tPTC5RiNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC5RiNCounts <- tPTC5RiNCounts[ptc5RiNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc5RiN <- CreateSeuratObject(f_tPTC5RiNCounts, meta.data = f_ptc5RiN@meta.data,
                                project = "ptc5rin")
# RE_ASSESS QC METRICS
f_ptc5RiN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc5RiN,file="processedData/filtered_PTC5_RightN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc5RiN_phase<-NormalizeData(f_ptc5RiN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc5RiN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc5RiN<-SCTransform(f_ptc5RiN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc5RiN<-f_ptc5RiN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc5RiN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc5RiN)<-"SCT_snn_res.0.6"
f_ptc5RiN<-f_ptc5RiN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc5RiN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc5RiN<-CellCycleScoring(f_ptc5RiN,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc5RiN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc5RiN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:13),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc5RiN_PCdata<-FetchData(f_ptc5RiN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc5RiN_umap_label <- FetchData(f_ptc5RiN,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:13), function(pc){
  ggplot(f_ptc5RiN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc5RiN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc5RiN[["pca"]], dims = 1:13, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc5RiN_markers<-FindAllMarkers(object = f_ptc5RiN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc5RiN_markers<-fptc5RiN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc5RiN_top10<-fptc5RiN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification

###########Thyroid 3 Right N##################
#load LN mets tissue dataset
ptc3RiN.data<-Read10X(data.dir = "scRNAseq Data/PTC3+RightN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc3RiN<-CreateSeuratObject(counts = ptc3RiN.data,project = "ptc3rin", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc3RiN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc3RiN$log10GenesPerUMI <- log10(ptc3RiN$nFeature_RNA) / log10(ptc3RiN$nCount_RNA)
# Compute percent mito ratio
ptc3RiN$mitoRatio <- PercentageFeatureSet(object = ptc3RiN, pattern = "^MT-")
ptc3RiN$mitoRatio <- ptc3RiN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc3RiN$riboRatio <- PercentageFeatureSet(object = ptc3RiN, pattern = "^RP[LS]")
ptc3RiN$riboRatio <- ptc3RiN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc3RiNMetadata<-ptc3RiN@meta.data
#add cell IDs to metadata
ptc3RiNMetadata$cells<-rownames(ptc3RiNMetadata)
#rename columns
ptc3RiNMetadata<-ptc3RiNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc3RiN@meta.data<-ptc3RiNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc3RiNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc3RiNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc3RiNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc3RiNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc3RiNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc3RiN <- subset(x = ptc3RiN, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC3RiNCounts<-GetAssayData(object = f_ptc3RiN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC3RiNNonzero <- tPTC3RiNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc3RiNkeep_genes <- Matrix::rowSums(tPTC3RiNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC3RiNCounts <- tPTC3RiNCounts[ptc3RiNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc3RiN <- CreateSeuratObject(f_tPTC3RiNCounts, meta.data = f_ptc3RiN@meta.data,
                                project = "ptc3rin")
# RE_ASSESS QC METRICS
f_ptc3RiN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc3RiN,file="processedData/filtered_PTC3_RightN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc3RiN_phase<-NormalizeData(f_ptc3RiN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc3RiN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc3RiN<-SCTransform(f_ptc3RiN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc3RiN<-f_ptc3RiN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc3RiN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc3RiN)<-"SCT_snn_res.0.6"
f_ptc3RiN<-f_ptc3RiN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc3RiN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc3RiN<-CellCycleScoring(f_ptc3RiN,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc3RiN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc3RiN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:18),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc3RiN_PCdata<-FetchData(f_ptc3RiN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc3RiN_umap_label <- FetchData(f_ptc3RiN,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:18), function(pc){
  ggplot(f_ptc3RiN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc3RiN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc3RiN[["pca"]], dims = 1:18, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc3RiN_markers<-FindAllMarkers(object = f_ptc3RiN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc3RiN_markers<-fptc3RiN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc3RiN_top10<-fptc3RiN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification

###########Thyroid 6 Right N####################
#load LN mets tissue dataset
ptc6RiN.data<-Read10X(data.dir = "scRNAseq Data/PTC6+RightN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc6RiN<-CreateSeuratObject(counts = ptc6RiN.data,project = "ptc6rin", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc6RiN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc6RiN$log10GenesPerUMI <- log10(ptc6RiN$nFeature_RNA) / log10(ptc6RiN$nCount_RNA)
# Compute percent mito ratio
ptc6RiN$mitoRatio <- PercentageFeatureSet(object = ptc6RiN, pattern = "^MT-")
ptc6RiN$mitoRatio <- ptc6RiN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc6RiN$riboRatio <- PercentageFeatureSet(object = ptc6RiN, pattern = "^RP[LS]")
ptc6RiN$riboRatio <- ptc6RiN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc6RiNMetadata<-ptc6RiN@meta.data
#add cell IDs to metadata
ptc6RiNMetadata$cells<-rownames(ptc6RiNMetadata)
#rename columns
ptc6RiNMetadata<-ptc6RiNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc6RiN@meta.data<-ptc6RiNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc6RiNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc6RiNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc6RiNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc6RiNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc6RiNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc6RiN <- subset(x = ptc6RiN, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC6RiNCounts<-GetAssayData(object = f_ptc6RiN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC6RiNNonzero <- tPTC6RiNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc6RiNkeep_genes <- Matrix::rowSums(tPTC6RiNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC6RiNCounts <- tPTC6RiNCounts[ptc6RiNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc6RiN <- CreateSeuratObject(f_tPTC6RiNCounts, meta.data = f_ptc6RiN@meta.data,
                                project = "ptc6rin")
# RE_ASSESS QC METRICS
f_ptc6RiN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc6RiN,file="processedData/filtered_PTC6_RightN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc6RiN_phase<-NormalizeData(f_ptc6RiN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc6RiN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc6RiN<-SCTransform(f_ptc6RiN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc6RiN<-f_ptc6RiN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc6RiN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc6RiN)<-"SCT_snn_res.0.6"
f_ptc6RiN<-f_ptc6RiN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc6RiN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc6RiN<-CellCycleScoring(f_ptc6RiN,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc6RiN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc6RiN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc6RiN_PCdata<-FetchData(f_ptc6RiN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc6RiN_umap_label <- FetchData(f_ptc6RiN,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc6RiN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc6RiN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc6RiN[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc6RiN_markers<-FindAllMarkers(object = f_ptc6RiN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc6RiN_markers<-fptc6RiN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc6RiN_top10<-fptc6RiN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification

###########Thyroid 7 Right N####################
#load LN mets tissue dataset
ptc7RiN.data<-Read10X(data.dir = "scRNAseq Data/PTC7+RightN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc7RiN<-CreateSeuratObject(counts = ptc7RiN.data,project = "ptc7rin", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc7RiN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc7RiN$log10GenesPerUMI <- log10(ptc7RiN$nFeature_RNA) / log10(ptc7RiN$nCount_RNA)
# Compute percent mito ratio
ptc7RiN$mitoRatio <- PercentageFeatureSet(object = ptc7RiN, pattern = "^MT-")
ptc7RiN$mitoRatio <- ptc7RiN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc7RiN$riboRatio <- PercentageFeatureSet(object = ptc7RiN, pattern = "^RP[LS]")
ptc7RiN$riboRatio <- ptc7RiN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc7RiNMetadata<-ptc7RiN@meta.data
#add cell IDs to metadata
ptc7RiNMetadata$cells<-rownames(ptc7RiNMetadata)
#rename columns
ptc7RiNMetadata<-ptc7RiNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc7RiN@meta.data<-ptc7RiNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc7RiNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc7RiNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc7RiNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc7RiNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc7RiNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc7RiN <- subset(x = ptc7RiN, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC7RiNCounts<-GetAssayData(object = f_ptc7RiN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC7RiNNonzero <- tPTC7RiNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc7RiNkeep_genes <- Matrix::rowSums(tPTC7RiNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC7RiNCounts <- tPTC7RiNCounts[ptc7RiNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc7RiN <- CreateSeuratObject(f_tPTC7RiNCounts, meta.data = f_ptc7RiN@meta.data,
                                project = "ptc7rin")
# RE_ASSESS QC METRICS
f_ptc7RiN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc7RiN,file="processedData/filtered_PTC7_RightN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc7RiN_phase<-NormalizeData(f_ptc7RiN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc7RiN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc7RiN<-SCTransform(f_ptc7RiN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc7RiN<-f_ptc7RiN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc7RiN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc7RiN)<-"SCT_snn_res.0.6"
f_ptc7RiN<-f_ptc7RiN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc7RiN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc7RiN<-CellCycleScoring(f_ptc7RiN,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc7RiN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc7RiN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc7RiN_PCdata<-FetchData(f_ptc7RiN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc7RiN_umap_label <- FetchData(f_ptc7RiN,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc7RiN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc6RiN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc7RiN[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc7RiN_markers<-FindAllMarkers(object = f_ptc7RiN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc7RiN_markers<-fptc7RiN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc7RiN_top10<-fptc7RiN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification


###########Thyroid 8 Tumor###############
#load tumor-pre dataset
ptc8T.data<-Read10X(data.dir = "scRNAseq Data/PTC8+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc8T<-CreateSeuratObject(counts = ptc8T.data,project = "ptc8t", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc8T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc8T$log10GenesPerUMI <- log10(ptc8T$nFeature_RNA) / log10(ptc8T$nCount_RNA)
# Compute percent mito ratio
ptc8T$mitoRatio <- PercentageFeatureSet(object = ptc8T, pattern = "^MT-")
ptc8T$mitoRatio <- ptc8T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc8T$riboRatio <- PercentageFeatureSet(object = ptc8T, pattern = "^RP[LS]")
ptc8T$riboRatio <- ptc8T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc8TMetadata<-ptc8T@meta.data
#add cell IDs to metadata
ptc8TMetadata$cells<-rownames(ptc8TMetadata)
#rename columns
ptc8TMetadata<-ptc8TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc8T@meta.data<-ptc8TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc8TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc8TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc8TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc8TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc8TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc8T <- subset(x = ptc8T, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC8tCounts<-GetAssayData(object = f_ptc8T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC8tNonzero <- tPTC8tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc8tkeep_genes <- Matrix::rowSums(tPTC8tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC8tCounts <- tPTC8tCounts[ptc8tkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc8T <- CreateSeuratObject(f_tPTC8tCounts, meta.data = f_ptc8T@meta.data,
                              project = "ptc8t")
# RE_ASSESS QC METRICS
f_ptc8T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc8T,file="processedData/filtered_PTC8_T.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc8T_phase<-NormalizeData(f_ptc8T)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc8T_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc8T<-SCTransform(f_ptc8T,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc8T<-f_ptc8T%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc8T@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc8T)<-"SCT_snn_res.0.6"
f_ptc8T<-f_ptc8T%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc8T,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc8T<-CellCycleScoring(f_ptc8T,
                          g2m.features = g2m_genes,
                          s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc8T,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc8T,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:13),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc8T_PCdata<-FetchData(f_ptc8T,vars = columns)
# Adding cluster label to center of cluster on UMAP
fp8t_umap_label <- FetchData(f_ptc8T,
                             vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:13), function(pc){
  ggplot(f_ptc8T_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fp8t_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc8T[["pca"]], dims = 1:13, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fp8t_markers<-FindAllMarkers(object = f_ptc8T,
                             only.pos = T,
                             logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fp8t_markers<-fp8t_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fp8t_top10<-fp8t_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
## Unsupervised Cluster Identification


###########Thyroid 9 Tumor###############
#load tumor-pre dataset
ptc9T.data<-Read10X(data.dir = "scRNAseq Data/PTC9+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc9T<-CreateSeuratObject(counts = ptc9T.data,project = "ptc9t", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc9T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc9T$log10GenesPerUMI <- log10(ptc9T$nFeature_RNA) / log10(ptc9T$nCount_RNA)
# Compute percent mito ratio
ptc9T$mitoRatio <- PercentageFeatureSet(object = ptc9T, pattern = "^MT-")
ptc9T$mitoRatio <- ptc9T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc9T$riboRatio <- PercentageFeatureSet(object = ptc9T, pattern = "^RP[LS]")
ptc9T$riboRatio <- ptc9T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc9TMetadata<-ptc9T@meta.data
#add cell IDs to metadata
ptc9TMetadata$cells<-rownames(ptc9TMetadata)
#rename columns
ptc9TMetadata<-ptc9TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc9T@meta.data<-ptc9TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc9TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc9TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc9TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc9TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc9TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc9T <- subset(x = ptc9T, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC9tCounts<-GetAssayData(object = f_ptc9T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC9tNonzero <- tPTC9tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc9tkeep_genes <- Matrix::rowSums(tPTC9tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC9tCounts <- tPTC9tCounts[ptc9tkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc9T <- CreateSeuratObject(f_tPTC9tCounts, meta.data = f_ptc9T@meta.data,
                              project = "ptc9t")
# RE_ASSESS QC METRICS
f_ptc9T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc9T,file="processedData/filtered_PTC9_T.RData")
#####
###########Thyroid 10 Tumor###############
#load tumor-pre dataset
ptc10T.data<-Read10X(data.dir = "scRNAseq Data/PTC10+T/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc10T<-CreateSeuratObject(counts = ptc10T.data,project = "ptc10t", 
                          min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc10T@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc10T$log10GenesPerUMI <- log10(ptc10T$nFeature_RNA) / log10(ptc10T$nCount_RNA)
# Compute percent mito ratio
ptc10T$mitoRatio <- PercentageFeatureSet(object = ptc10T, pattern = "^MT-")
ptc10T$mitoRatio <- ptc10T@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc10T$riboRatio <- PercentageFeatureSet(object = ptc10T, pattern = "^RP[LS]")
ptc10T$riboRatio <- ptc10T@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc10TMetadata<-ptc10T@meta.data
#add cell IDs to metadata
ptc10TMetadata$cells<-rownames(ptc10TMetadata)
#rename columns
ptc10TMetadata<-ptc10TMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc10T@meta.data<-ptc10TMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc10TMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc10TMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc10TMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc10TMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc10TMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc10T <- subset(x = ptc10T, 
                  subset= (nUMI >= 500) & 
                    (nGene >= 250) & 
                    (log10GenesPerUMI > 0.80) & 
                    (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC10tCounts<-GetAssayData(object = f_ptc10T,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC10tNonzero <- tPTC10tCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc10tkeep_genes <- Matrix::rowSums(tPTC10tNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC10tCounts <- tPTC10tCounts[ptc10tkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc10T <- CreateSeuratObject(f_tPTC10tCounts, meta.data = f_ptc10T@meta.data,
                              project = "ptc10t")
# RE_ASSESS QC METRICS
f_ptc10T@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc10T,file="processedData/filtered_PTC10_T.RData")
#####
###########Thyroid 10 Right N####################
#load LN mets tissue dataset
ptc10RiN.data<-Read10X(data.dir = "scRNAseq Data/PTC10+RightN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc10RiN<-CreateSeuratObject(counts = ptc10RiN.data,project = "ptc10rin", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc10RiN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc10RiN$log10GenesPerUMI <- log10(ptc10RiN$nFeature_RNA) / log10(ptc10RiN$nCount_RNA)
# Compute percent mito ratio
ptc10RiN$mitoRatio <- PercentageFeatureSet(object = ptc10RiN, pattern = "^MT-")
ptc10RiN$mitoRatio <- ptc10RiN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc10RiN$riboRatio <- PercentageFeatureSet(object = ptc10RiN, pattern = "^RP[LS]")
ptc10RiN$riboRatio <- ptc10RiN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc10RiNMetadata<-ptc10RiN@meta.data
#add cell IDs to metadata
ptc10RiNMetadata$cells<-rownames(ptc10RiNMetadata)
#rename columns
ptc10RiNMetadata<-ptc10RiNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc10RiN@meta.data<-ptc10RiNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc10RiNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc10RiNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc10RiNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc10RiNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc10RiNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc10RiN <- subset(x = ptc10RiN, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC10RiNCounts<-GetAssayData(object = f_ptc10RiN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC10RiNNonzero <- tPTC10RiNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc10RiNkeep_genes <- Matrix::rowSums(tPTC10RiNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC10RiNCounts <- tPTC10RiNCounts[ptc10RiNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc10RiN <- CreateSeuratObject(f_tPTC10RiNCounts, meta.data = f_ptc10RiN@meta.data,
                                project = "ptc10rin")
# RE_ASSESS QC METRICS
f_ptc10RiN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc10RiN,file="processedData/filtered_PTC10_RightN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc10RiN_phase<-NormalizeData(f_ptc10RiN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc10RiN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc10RiN<-SCTransform(f_ptc10RiN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc10RiN<-f_ptc10RiN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc10RiN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc10RiN)<-"SCT_snn_res.0.6"
f_ptc10RiN<-f_ptc10RiN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc10RiN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc10RiN<-CellCycleScoring(f_ptc10RiN,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc10RiN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc10RiN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc10RiN_PCdata<-FetchData(f_ptc10RiN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc10RiN_umap_label <- FetchData(f_ptc10RiN,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc10RiN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc10RiN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc10RiN[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc10RiN_markers<-FindAllMarkers(object = f_ptc10RiN,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc10RiN_markers<-fptc10RiN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc10RiN_top10<-fptc10RiN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification
#####
###########Thyroid 11 Right N####################
#load LN mets tissue dataset
ptc11RiN.data<-Read10X(data.dir = "scRNAseq Data/PTC11+RightN/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc11RiN<-CreateSeuratObject(counts = ptc11RiN.data,project = "ptc11rin", 
                             min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc11RiN@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc11RiN$log10GenesPerUMI <- log10(ptc11RiN$nFeature_RNA) / log10(ptc11RiN$nCount_RNA)
# Compute percent mito ratio
ptc11RiN$mitoRatio <- PercentageFeatureSet(object = ptc11RiN, pattern = "^MT-")
ptc11RiN$mitoRatio <- ptc11RiN@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc11RiN$riboRatio <- PercentageFeatureSet(object = ptc11RiN, pattern = "^RP[LS]")
ptc11RiN$riboRatio <- ptc11RiN@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc11RiNMetadata<-ptc11RiN@meta.data
#add cell IDs to metadata
ptc11RiNMetadata$cells<-rownames(ptc11RiNMetadata)
#rename columns
ptc11RiNMetadata<-ptc11RiNMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc11RiN@meta.data<-ptc11RiNMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc11RiNMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc11RiNMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc11RiNMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc11RiNMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc11RiNMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc11RiN <- subset(x = ptc11RiN, 
                     subset= (nUMI >= 500) & 
                       (nGene >= 250) & 
                       (log10GenesPerUMI > 0.80) & 
                       (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC11RiNCounts<-GetAssayData(object = f_ptc11RiN,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC11RiNNonzero <- tPTC11RiNCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc11RiNkeep_genes <- Matrix::rowSums(tPTC11RiNNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC11RiNCounts <- tPTC11RiNCounts[ptc11RiNkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc11RiN <- CreateSeuratObject(f_tPTC11RiNCounts, meta.data = f_ptc11RiN@meta.data,
                                 project = "ptc11rin")
# RE_ASSESS QC METRICS
f_ptc11RiN@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc11RiN,file="processedData/filtered_PTC11_RightN.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc11RiN_phase<-NormalizeData(f_ptc11RiN)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc11RiN_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc11RiN<-SCTransform(f_ptc11RiN,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc11RiN<-f_ptc11RiN%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc11RiN@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc11RiN)<-"SCT_snn_res.0.6"
f_ptc11RiN<-f_ptc11RiN%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc11RiN,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc11RiN<-CellCycleScoring(f_ptc11RiN,
                             g2m.features = g2m_genes,
                             s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc11RiN,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc11RiN,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc11RiN_PCdata<-FetchData(f_ptc11RiN,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc11RiN_umap_label <- FetchData(f_ptc11RiN,
                                  vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc11RiN_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc11RiN_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc11RiN[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc11RiN_markers<-FindAllMarkers(object = f_ptc11RiN,
                                  only.pos = T,
                                  logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc11RiN_markers<-fptc11RiN_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc11RiN_top10<-fptc11RiN_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification
#####
###########Thyroid 4 Subcutaneous####################
#load LN mets tissue dataset
ptc4SC.data<-Read10X(data.dir = "scRNAseq Data/PTC4+SC/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc4SC<-CreateSeuratObject(counts = ptc4SC.data,project = "ptc4sc", 
                             min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc4SC@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc4SC$log10GenesPerUMI <- log10(ptc4SC$nFeature_RNA) / log10(ptc4SC$nCount_RNA)
# Compute percent mito ratio
ptc4SC$mitoRatio <- PercentageFeatureSet(object = ptc4SC, pattern = "^MT-")
ptc4SC$mitoRatio <- ptc4SC@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc4SC$riboRatio <- PercentageFeatureSet(object = ptc4SC, pattern = "^RP[LS]")
ptc4SC$riboRatio <- ptc4SC@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc4SCMetadata<-ptc4SC@meta.data
#add cell IDs to metadata
ptc4SCMetadata$cells<-rownames(ptc4SCMetadata)
#rename columns
ptc4SCMetadata<-ptc4SCMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc4SC@meta.data<-ptc4SCMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc4SCMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc4SCMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc4SCMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc4SCMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc4SCMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc4SC <- subset(x = ptc4SC, 
                     subset= (nUMI >= 500) & 
                       (nGene >= 250) & 
                       (log10GenesPerUMI > 0.80) & 
                       (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC4SCCounts<-GetAssayData(object = f_ptc4SC,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC4SCNonzero <- tPTC4SCCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc4SCkeep_genes <- Matrix::rowSums(tPTC4SCNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC4SCCounts <- tPTC4SCCounts[ptc4SCkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc4SC <- CreateSeuratObject(f_tPTC4SCCounts, meta.data = f_ptc4SC@meta.data,
                                 project = "ptc4sc")
# RE_ASSESS QC METRICS
f_ptc4SC@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc4SC,file="processedData/filtered_PTC4_Subcutaneous.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc4SC_phase<-NormalizeData(f_ptc4SC)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc4SC_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc4SC<-SCTransform(f_ptc4SC,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc4SC<-f_ptc4SC%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc4SC@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc4SC)<-"SCT_snn_res.0.6"
f_ptc4SC<-f_ptc4SC%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc4SC,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc4SC<-CellCycleScoring(f_ptc4SC,
                             g2m.features = g2m_genes,
                             s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc4SC,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc4SC,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc4SC_PCdata<-FetchData(f_ptc4SC,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc4SC_umap_label <- FetchData(f_ptc4SC,
                                  vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc4SC_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc4SC_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc4SC[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc4SC_markers<-FindAllMarkers(object = f_ptc4SC,
                                  only.pos = T,
                                  logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc4SC_markers<-fptc4SC_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc4SC_top10<-fptc4SC_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification
#####
###########Thyroid 11 Subcutaneous####################
#load LN mets tissue dataset
ptc11SC.data<-Read10X(data.dir = "scRNAseq Data/PTC11+SC/")
#initialize the seurat object with the raw (non-normalized data)
#remove the genes detected in less than 3 cells and cell samples expressing less than 200 genes
ptc11SC<-CreateSeuratObject(counts = ptc11SC.data,project = "ptc11sc", 
                            min.cells = 3,min.features = 200)

## QUALITY CONTROL
head(ptc11SC@meta.data)
# Add number of genes per UMI for each cell to metadata
# more genes detected per UMI, more complex our data
ptc11SC$log10GenesPerUMI <- log10(ptc11SC$nFeature_RNA) / log10(ptc11SC$nCount_RNA)
# Compute percent mito ratio
ptc11SC$mitoRatio <- PercentageFeatureSet(object = ptc11SC, pattern = "^MT-")
ptc11SC$mitoRatio <- ptc11SC@meta.data$mitoRatio / 100
# Compute percent ribo ratio
ptc11SC$riboRatio <- PercentageFeatureSet(object = ptc11SC, pattern = "^RP[LS]")
ptc11SC$riboRatio <- ptc11SC@meta.data$riboRatio / 100

# GENERATING QUALITY METRICS
# Additional metadata column
ptc11SCMetadata<-ptc11SC@meta.data
#add cell IDs to metadata
ptc11SCMetadata$cells<-rownames(ptc11SCMetadata)
#rename columns
ptc11SCMetadata<-ptc11SCMetadata%>%
  dplyr::rename(seq_folder=orig.ident,
                nUMI=nCount_RNA,
                nGene=nFeature_RNA)
#add metatdata back to Seurat object
ptc11SC@meta.data<-ptc11SCMetadata

# ASSESSING THE QUALITY METRICS
# Cell count: check the objects number of ptc1T
# UMI counts (transcripts) per cell
# The UMI counts per cell should generally be above 500
ptc11SCMetadata%>%
  ggplot(aes(x=nUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 500)
# Genes detected per cell
# For high quality data, the proportional histogram should 
# contain a single large peak that represents cells that were encapsulated.
ptc11SCMetadata%>%
  ggplot(aes(x=nGene))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  ylab("Cell density")+
  geom_vline(xintercept = 300)
# Complexity
# Generally, we expect the novelty score to be above 0.80 for good quality cells.
ptc11SCMetadata%>%
  ggplot(aes(x=log10GenesPerUMI))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.8)
# Mitochondrial counts ratio
# poor quality samples (dead or dying cells) for mitochondrial counts 
# as cells which surpass the 0.2 mitochondrial ratio mark
ptc11SCMetadata%>%
  ggplot(aes(x=mitoRatio))+
  geom_density()+
  scale_x_log10()+
  theme_classic()+
  geom_vline(xintercept = 0.2)

# Joint filtering effects
# Visualize the correlation between genes detected and number of UMIs and determine 
# whether strong presence of cells with low numbers of genes/UMIs
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs 
# (upper right quadrant of the plot).
# Mitochondrial read fractions are only high in particularly low count cells 
# with few detected genes (darker colored data points)
ptc11SCMetadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# FILTERING
# Cell-level filtering
# nUMI>500, nGene>250, log10GenesPerUMI>0.8, mitoRatio<0.2
# Filter out low quality cells using selected thresholds - these will change with experiment
f_ptc11SC <- subset(x = ptc11SC, 
                    subset= (nUMI >= 500) & 
                      (nGene >= 250) & 
                      (log10GenesPerUMI > 0.80) & 
                      (mitoRatio < 0.20))
# Gene-level filtering
# remove genes with zero counts
tPTC11SCCounts<-GetAssayData(object = f_ptc11SC,slot = "counts")
# Output a logical matrix specifying for each gene on whether or 
# not there are more than zero counts per cell
tPTC11SCNonzero <- tPTC11SCCounts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
ptc11SCkeep_genes <- Matrix::rowSums(tPTC11SCNonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
f_tPTC11SCCounts <- tPTC11SCCounts[ptc11SCkeep_genes, ]
# Reassign to filtered Seurat object
f_ptc11SC <- CreateSeuratObject(f_tPTC11SCCounts, meta.data = f_ptc11SC@meta.data,
                                project = "ptc11sc")
# RE_ASSESS QC METRICS
f_ptc11SC@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# SAVE
save(f_ptc11SC,file="processedData/filtered_PTC11_Subcutaneous.RData")

## SCT NORMALIZATION

# Normalize the counts
# load cycle.rda
# Score cells for cell cycle
# identify the most variable genes
# scale the counts, perform PCA, plot the PCA colored by cell cycle phase
f_ptc11SC_phase<-NormalizeData(f_ptc11SC)%>%
  CellCycleScoring(g2m.features = g2m_genes,
                   s.features = s_genes)%>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000,
                       verbose = T)%>%
  ScaleData()%>%
  RunPCA()
# Visualization
# No significant difference in pattern between different cell phases
DimPlot(f_ptc11SC_phase,reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

# SCTransform to Normalization and regressing out unwanted variation
f_ptc11SC<-SCTransform(f_ptc11SC,vars.to.regress = c("mitoRatio"))

## CLUSTERING
# Find Neighbors
#Determine the K-nearest neighbor graph
f_ptc11SC<-f_ptc11SC%>%
  RunPCA()%>%
  FindNeighbors(dims = 1:40)%>%
  FindClusters(resolution = c(0.6,0.8,1.0,1.4))
# check the result
View(f_ptc11SC@meta.data)

# Visualize clusters of cells
# Explore resolutions
#Assign identity of clusters
#usually start with 0.6-0.8
Idents(object = f_ptc11SC)<-"SCT_snn_res.0.6"
f_ptc11SC<-f_ptc11SC%>%
  RunUMAP(reduction = "pca",
          dims = 1:40)
DimPlot(f_ptc11SC,
        reduction = "umap",
        label = T,
        label.size = 6,
        repel = T)
# incorporate the cell cycle score
f_ptc11SC<-CellCycleScoring(f_ptc11SC,
                            g2m.features = g2m_genes,
                            s.features = s_genes)

## CLUSTER IDENTIFICATION

# Exploration of quality control metrics
# Segregation of clusters by cell cycle phase
DimPlot(f_ptc11SC,reduction = "pca",
        label = T,
        split.by = "Phase")
#The order argument will plot the positive cells above the negative cells, 
#while the min.cutoff argument will determine the threshold for shading. 
#A min.cutoff of q10 translates to the 10% of cells with the lowest expression of the gene 
#will not exhibit any purple shading (completely gray).
metrics<-c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
FeaturePlot(f_ptc11SC,
            reduction = "umap",
            features = metrics,
            pt.size = 0.4,
            order = T,
            min.cutoff = 'q10',
            label = T)

#Exploration of the PCs driving the different clusters
#Defining the information in the seurat object of interest
columns<-c(paste0("PC_",1:17),
           "ident",
           "UMAP_1","UMAP_2")
#Extracting this data from the seurat object
f_ptc11SC_PCdata<-FetchData(f_ptc11SC,vars = columns)
# Adding cluster label to center of cluster on UMAP
fptc11SC_umap_label <- FetchData(f_ptc11SC,
                                 vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))
# Plotting a UMAP plot for each of the PCs
# Try to identify the similarity and uniqueness between cluster
map(paste0("PC_", 1:17), function(pc){
  ggplot(f_ptc11SC_PCdata, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=fptc11SC_umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results
print(f_ptc11SC[["pca"]], dims = 1:17, nfeatures = 10)
## MARKER IDENTIFICATION
# To identify the gene markers for each cluster and identify cell types with them
# Identification of all markers for each cluster
fptc11SC_markers<-FindAllMarkers(object = f_ptc11SC,
                                 only.pos = T,
                                 logfc.threshold = 0.25)
# Add annotation
annotations<-read.csv("annotation.csv")
# Combine marker with gene description
fptc11SC_markers<-fptc11SC_markers%>%
  left_join(y = unique(annotations[,c("gene_name","description")]),
            by = c("gene" = "gene_name"))
# Extract top 10 markers per cluster
fptc11SC_top10<-fptc11SC_markers%>%
  group_by(cluster)%>%
  top_n(n = 10, wt = avg_log2FC)
# Unsupervised Cluster Identification
#####
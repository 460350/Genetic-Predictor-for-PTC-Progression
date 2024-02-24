library(Seurat)
library(SCEVAN)
library(tidyverse)

setwd('C:/Users/Lenovo/Documents/R Data/PTC_Pair')

####ptc11sc####
load(file = "combined_ptcnonImmune_cell.RData")
Idents(object = fh_NonImmune)<-"orig.ident"
ptc11sc<-subset(fh_NonImmune,idents= 'ptc11sc')
rm(fh_NonImmune)

ptc11sc_mtx<-ptc11sc@assays$RNA@counts

# Time consuming steps: [10] Adjust baseline
# results can be slightly different when set different par_cores
ptc11sc.results <- SCEVAN::pipelineCNA(ptc11sc_mtx, sample = "ptc11sc",
                                     par_cores = 10,
                                     SUBCLONES = F, plotTree = F)

View(ptc11sc.results)

ptc11sc<-AddMetaData(ptc11sc,metadata = ptc11sc.results,
                   col.name = c("SCEVANpred","SCEVANconfi"))

####SAVE####
save(ptc11sc,file="Split Sample/ptc11sc.RData")

rm(ptc11sc)
#####


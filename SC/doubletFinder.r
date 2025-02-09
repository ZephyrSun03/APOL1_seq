
## R/4.0.3

Args <- commandArgs(trailingOnly = T)

library(Seurat)
library(DoubletFinder)
library(RColorBrewer)
library(ggplot2)

ObjF <- Args[1]

Colors <- c(brewer.pal(8, "Set2")[c(-6, -8)], brewer.pal(12, "Set3")[c(-2, -3)], brewer.pal(9, "Set1")[c(-6, -9)], brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))

Obj.all <- readRDS(file = ObjF)

CName <- colnames(Obj.all@meta.data)

Sams <- unique(Obj.all$orig.ident)

ObjList <- list()

for(ss in Sams){

Obj <- subset(Obj.all, cells = colnames(Obj.all)[Obj.all$orig.ident == ss])

pdf(paste0(ss, ".pK.distribution.pdf"))

	sweep.res.list <- paramSweep_v3(Obj, PCs = 1:10, sct = FALSE)
	
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	
	bcmvn <- find.pK(sweep.stats)
	
	pKNum <- as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))

dev.off()

#stop()

## Homotypic Doublet Proportion Estimate ---------------------------------------------------------------------

homotypic.prop <- modelHomotypic(Obj@meta.data$category_anno_meta)    ### ex: annotations <- Obj@meta.data$ClusteringResults

nExp_poi <- round(homotypic.prop*nrow(Obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Run DoubletFinder with varying classification stringencies ------------------------------------------------

Obj <- doubletFinder_v3(Obj, PCs = 1:10, pN = 0.25, pK = pKNum, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

#saveRDS(Obj, file = "Obj.DB.rds")

Obj <- doubletFinder_v3(Obj, PCs = 1:10, pN = 0.25, pK = pKNum, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_", pKNum, "_", nExp_poi), sct = FALSE)

ObjList[[ss]] <- Obj

}

#saveRDS(ObjList, file = "ObjDBList.rds")

DBTab <- Reduce(rbind, lapply(ObjList, function(x) {Tab <- x@meta.data; colnames(Tab) <- c(CName, "pAnn", "DF1", "DF2"); Tab}))

write.table(DBTab, file = "DBTab.xls", sep = "\t", quote = F, col.names = NA)


Obj.all$DF1 <- DBTab[colnames(Obj.all), "DF1"]

Obj.all$DF2 <- DBTab[colnames(Obj.all), "DF2"]

pdf(paste0("Doublet1", "_TSNEPlot.pdf"), width = 10)

        DimPlot(Obj.all, reduction = "umap", label = F, repel = T, label.size = 4, pt.size = 0.3, group.by = "DF1", cols = Colors) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()


pdf(paste0("Doublet2", "_TSNEPlot.pdf"), width = 10)

        DimPlot(Obj.all, reduction = "umap", label = F, repel = T, label.size = 4, pt.size = 0.3, group.by = "DF2", cols = Colors) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()


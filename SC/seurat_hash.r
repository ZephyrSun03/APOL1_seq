
## R/4.0.3

Args <- commandArgs(trailingOnly = T)

library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

Colors <- c(brewer.pal(8, "Set2")[c(-6, -8)], brewer.pal(9, "Set1")[c(-1, -6, -9)], brewer.pal(12, "Set3")[c(-2)], brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))


CellF <- Args[1]
ObjS <- Args[2]
QuantN <- as.numeric(Args[3])
KeepS <- Args[4]

Mat <- Read10X(data.dir = CellF)
Obj <- readRDS(file = ObjS)


Mat <- Mat[["Antibody Capture"]]

colnames(Mat) <- paste0("APOL1_", colnames(Mat))

apply(Mat, 1, function(x) quantile(x, probs = seq(0, 1, by = 0.1)))

apply(Mat, 1, function(x) quantile(x[x > 0], probs = seq(0, 1, by = 0.1)))

rowSums(Mat)


Obj[["HTO"]] <- CreateAssayObject(counts = Mat[, colnames(Obj)])

Obj <- NormalizeData(Obj, assay = "HTO", normalization.method = "CLR")

Mat.norm <- GetAssayData(Obj, slot = "data", assay = "HTO")

apply(Mat.norm, 1, function(x) quantile(x, probs = seq(0, 1, by = 0.1)))

apply(Mat.norm, 1, function(x) quantile(x[x > 0], probs = seq(0, 1, by = 0.1)))

rowSums(Mat.norm)



Obj <- HTODemux(Obj, assay = "HTO", positive.quantile = QuantN)

#Obj <- HTODemux(Obj, assay = "HTO", positive.quantile = 0.95)

table(Obj$HTO_classification.global)

table(Obj$HTO_maxID)

table(Obj@meta.data[, c("HTO_maxID", "HTO_classification.global")])

table(Obj@meta.data[, c("category_anno_meta", "HTO_classification.global")])

#write.table(table(Obj@meta.data[, c("category_anno_meta", "HTO_maxID")]), file = paste0(ObjS, "_Sample_Ann.xls"), sep = "\t", quote = F, col.names = NA)

# Group cells based on the max HTO signal

Obj$orig.ident <- Obj$HTO_maxID 

Idents(Obj) <- "HTO_maxID"

pdf("HTO_maxID.data.pdf", height = 14)

	RidgePlot(Obj, assay = "HTO", slot = "data", features = rownames(Obj[["HTO"]]), ncol = 1)

dev.off()

pdf("HTO_maxID.counts.pdf", height = 14)

        RidgePlot(Obj, assay = "HTO", slot = "counts", features = rownames(Obj[["HTO"]]), ncol = 1)

dev.off()

pdf("HTO_maxID.logcounts.pdf", height = 14)

        RidgePlot(Obj, assay = "HTO", slot = "counts", features = rownames(Obj[["HTO"]]), ncol = 1, log = T)

dev.off()


Idents(Obj) <- "HTO_classification.global"

pdf("HTO_classification.global.pdf")

	VlnPlot(Obj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

dev.off()


pdf("HTO_heat.pdf", width = 12, height = 6)

	HTOHeatmap(Obj, assay = "HTO", ncells = 5000, raster = FALSE)

dev.off()


pdf(paste0(ObjS, "_HTO_classification.global", "_TSNEPlot_Ann.pdf"), width = 10)

        DimPlot(Obj, reduction = "umap", label = TRUE, repel = T, label.size = 4, pt.size = 0.3, group.by = "HTO_classification.global", cols = Colors) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()


Stat <- table(Obj@meta.data[, c("category_anno_meta", "HTO_classification.global")])

Stat <- cbind(Stat, DoubleRate = Stat[, "Doublet"]/rowSums(Stat))

write.table(Stat, file = "Stat_DoubleRate.xls", sep = "\t", quote = F, col.names = NA)


Obj <- subset(Obj, subset = HTO_classification.global %in% unlist(strsplit(KeepS, split = ";")))

write.table(table(Obj@meta.data[, c("category_anno_meta", "HTO_maxID")]), file = paste0(ObjS, "_Sample_Ann.xls"), sep = "\t", quote = F, col.names = NA)

pdf(paste0(ObjS, "_HTO_maxID", "_TSNEPlot_Ann.pdf"), width = 10)

        DimPlot(Obj, reduction = "umap", label = TRUE, repel = T, label.size = 4, pt.size = 0.3, group.by = "HTO_maxID", cols = Colors) + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()


Colors <- Colors[1:length(levels(Obj$category_anno_meta))]

names(Colors) <- levels(Obj$category_anno_meta)

plist <- list()

for(ss in unique(Obj$HTO_maxID)){

        ObjTmp <- subset(Obj, cells = names(Obj$HTO_maxID[Obj$HTO_maxID== ss]))

        pp  <- DimPlot(ObjTmp, reduction = "umap", label = TRUE, repel = T, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(ObjTmp$category_anno_meta)]) + NoLegend() + ggtitle(ss)

        plist[[ss]] <- pp

}

ggsave(file = paste0(ObjS, "_TSNEPlot_HashTag.pdf"), arrangeGrob(grobs = plist, ncol = 2), width = 12, height = 4 * (as.integer(length(plist)/2)+1), limitsize = FALSE)

Obj$orig.ident <- Obj$HTO_maxID

Idents(Obj) <- Obj$category_anno_meta


saveRDS(Obj, file = paste0(ObjS, ".hto.rds"))

## First, we will remove negative cells from the object
#
#Obj.subset <- subset(Obj, idents = "Negative", invert = TRUE)
#
## Calculate a tSNE embedding of the HTO data
#
#DefaultAssay(Obj.subset) <- "HTO"
#
#Obj.subset <- ScaleData(Obj.subset, features = rownames(Obj.subset),verbose = FALSE)
#
#Obj.subset <- RunPCA(Obj.subset, features = rownames(Obj.subset), approx = FALSE)
#
#Obj.subset <- RunTSNE(Obj.subset, dims = 1:6, perplexity = 100)
#
#pdf("HTO_classification.global.UMAP.pdf")
#
#	DimPlot(Obj.subset)
#
#dev.off()


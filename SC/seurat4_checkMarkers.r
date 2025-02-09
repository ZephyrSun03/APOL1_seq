
#! R/4.0.3

library(Seurat)
library(RColorBrewer)
library(lattice)
library(ggplot2)
library(gridExtra)
library(cowplot)

source("/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/customize_Seurat_FeaturePlot.r")

Args <- commandArgs(trailingOnly = TRUE)

ObjF <- Args[1]
KeyS <- Args[2]
AnnF <- Args[3]
MarkerF <- Args[4]
AssayS <- Args[5]
OrderS <- Args[6]
SepS <- Args[7]

Obj <- readRDS(ObjF)

category_anno <- read.table(file = AnnF, sep = "\t")

category_anno <- data.frame(category_anno[category_anno[[1]] %in% colnames(Obj),])

category_anno_meta <- data.frame(category_anno_meta = category_anno[, 3])

row.names(category_anno_meta) <- category_anno[[1]]

Obj <- subset(Obj, cells = rownames(category_anno_meta))

Obj <- AddMetaData(object = Obj, metadata = category_anno_meta, col.name = "category_anno_meta")

write.table(table(data.frame(Ann = Obj$category_anno_meta, Sample = Obj$orig.ident)), file = paste0(KeyS, "Sample_Ann.xls"), sep = "\t", quote = FALSE, col.names = NA)

#Expr <- GetAssayData(Obj, slot = "data")

#Colors <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))

Colors <- c(brewer.pal(8, "Set2")[c(-1, -6, -8)], brewer.pal(9, "Set1")[c(-1, -6, -9)], brewer.pal(12, "Set3")[c(-2)], brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))

Colors2 <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", header = F, row.names = 1, as.is = T)


pdf(paste0("14.TSNEPlot_Ann", KeyS, ".pdf"), width = 11)

	DimPlot(Obj, reduction = "umap", label = TRUE, pt.size = 0.3, label.size = 6, group.by = "category_anno_meta", cols = Colors[1:length(unique(category_anno[[3]]))]) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.position = "none")

dev.off()

#pdf("15.TSNEPlot_Ann_split.pdf", width = 16)
#
#	DimPlot(Obj, reduction = "umap", label = TRUE, pt.size = 0.3, group.by = "category_anno_meta", cols = brewer.pal(n = length(unique(category_anno[[3]])), name = "Set3"), split.by = "orig.ident")
#
#dev.off()


DefaultAssay(Obj) <- AssayS

if(OrderS != "Null"){

	Obj <- subset(Obj, idents = unlist(strsplit(OrderS, split = ";")))

	Idents(Obj) <- factor(as.character(Idents(Obj)), levels = rev(unlist(strsplit(OrderS, split = ";"))))

}else{

        Idents(Obj) <- Idents(Obj)
}


if(SepS == "Sep"){

        Idents(Obj) <- paste(Idents(Obj), Obj$orig.ident, sep = "_")

}


markers_c <- read.table(file = MarkerF, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

markers_cc <- intersect(markers_c$Gene, rownames(Obj))

#markers_cc <- factor(markers_cc, levels = markers_cc)

#markers_cc <- markers_cc[1:6]

        #plist <- FeaturePlot(object = Obj, features = markers_cc, cols = brewer.pal(9, "YlOrRd")[c(2, 9)], reduction = "umap", pt.size = 0.1, do.return = T)

	#plist <- customize_Seurat_FeaturePlot(p = plist, expression.threshold = 0)

	#ggsave(plist, file = paste0("11.Feature_plot_", Args[2], ".pdf"), width = 22, height = 6 * (as.integer(length(markers_cc)/3)+1), limitsize = FALSE)

        plist <- FeaturePlot(object = Obj, features = (markers_cc), cols = brewer.pal(9, "YlOrRd")[c(2, 9)], reduction = "umap", pt.size = 0.1, combine = FALSE)

	plist <- lapply(plist, function(x) x + theme(plot.title = element_text(size = 20)))

	#plist <- lapply(plist, function(x) customize_Seurat_FeaturePlot(p = x, expression.threshold = 0))

        ggsave(file = paste0("11.Feature_plot_", KeyS, ".pdf"), arrangeGrob(grobs = plist, ncol = 3), width = 22, height = 6 * (as.integer(length(markers_cc)/3)+1), limitsize = FALSE)


        plist <- VlnPlot(object = Obj, features = markers_cc, combine = FALSE, pt.size = 0, cols = Colors)

	plist <- lapply(plist, function(x) x + geom_boxplot(outlier.shape = NA) + theme(legend.position = "none", axis.title.x = element_blank()))

        ggsave(file = paste0("12.Vln_plot_", KeyS, ".pdf"), arrangeGrob(grobs = plist, ncol = 1), width = 0.6 * length(unique(Obj$category_anno_meta)), height = 2.5 * length(markers_cc), limitsize = FALSE)

markers_cc <- factor((markers_cc), levels = rev(markers_cc))

pdf(paste0("13.Dot_plot_", KeyS, ".pdf"), width = 0.3 *length(markers_cc) + 6, height = 0.25 * length(unique(Obj$category_anno_meta)) + 3)

	DotPlot(Obj, assay = "RNA", features = markers_cc, cols = Colors2[c("HeatPurple2", "HeatPurple8"), ], dot.scale = 8, scale = TRUE) + RotatedAxis() + theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size = 13), axis.title.y = element_blank(), axis.title.x = element_blank())

dev.off()

all.genes <- rownames(Obj)

Obj <- ScaleData(Obj, features = all.genes)

pp <- DoHeatmap(Obj, features = as.character(markers_cc), slot = "scale.data", size = 3.5, raster = FALSE, group.colors = Colors) + theme(axis.text.y = element_text(size = 10)) + scale_fill_gradientn(colors = brewer.pal(n = 11, name = "RdBu")[c(9, 6, 3)])

ggsave(paste0(KeyS, "_heatmap.pdf"), pp, height = 7, width = 12)


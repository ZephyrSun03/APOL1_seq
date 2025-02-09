
#! R/4.0.3

Args <- commandArgs(trailingOnly = TRUE)

library(Seurat)
library(RColorBrewer)
library(lattice)
library(ggplot2)
library(gridExtra)
library(dplyr)

ObjF <- Args[1]
MarkerL <- Args[2]
AnnF <- Args[3]
SamInfoF <- Args[4]
CellOrderS <- Args[5]

#Colors <- c(brewer.pal(8, "Set2")[c(-6, -8)], brewer.pal(12, "Set3")[c(-2, -3)], brewer.pal(9, "Set1")[c(-6, -9)], brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))

Colors <- c(brewer.pal(8, "Set2")[c(-6, -8)], brewer.pal(9, "Set1")[c(-1, -6, -9)], brewer.pal(12, "Set3")[c(-2)], brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"))

Obj <- readRDS(file = ObjF)

Obj$log10_nGene <- log10(Obj$nFeature_RNA + 1)

Obj$log10_nUMI <- log10(Obj$nCount_RNA + 1)

category_anno <- read.table(file = AnnF, sep = "\t")

category_anno <- data.frame(category_anno[category_anno[[1]] %in% colnames(Obj),])

category_anno_meta <- data.frame(category_anno_meta = category_anno[, 3])

row.names(category_anno_meta) <- category_anno[[1]]

Obj <- subset(Obj, cells = rownames(category_anno_meta))

Obj <- AddMetaData(object = Obj, metadata = category_anno_meta, col.name = "category_anno_meta")

if(CellOrderS !="Null"){

        Obj$category_anno_meta <- factor(as.character(Obj$category_anno_meta), levels = unlist(strsplit(CellOrderS, split = ";")))

}else{

	Obj$category_anno_meta <- factor(Obj$category_anno_meta)

}

SamInfo <- read.delim(file = SamInfoF, row.names = 1, stringsAsFactors = FALSE, as.is = T)

Obj <- subset(Obj, subset = orig.ident %in% rownames(SamInfo))

Obj$SamInfo <- factor(SamInfo[Obj$orig.ident, ], levels = SamInfo$Treatment[!duplicated(SamInfo$Treatment)])


SamAnn <- table(data.frame(Ann = Obj$category_anno_meta, Sample = Obj$orig.ident))

write.table(SamAnn, file = paste0(AnnF, "_Sample_Ann.xls"), sep = "\t", quote = FALSE, col.names = NA)

Idents(Obj) <- Obj$category_anno_meta

saveRDS(Obj, file = paste0(AnnF, "_Obj_post.rds"))

DefaultAssay(Obj) <- "RNA"

if(MarkerL == "Yes"){

        Obj.markers <- FindAllMarkers(Obj, only.pos = TRUE, min.pct = 0, logfc.threshold = 0, return.thresh = 1, test.use = "wilcox")

        write.table(Obj.markers, file = paste0(AnnF, "_Markers.lst.ann"), sep = "\t", quote = F, col.names = NA)

        print(paste(dim(Obj.markers)[1], "markers found!"));

	#Idents(Obj) <- factor(as.character(Idents(Obj)), c("Trem2_hi_Mac", "Resident_Mac", "Inflammatory_Mac", "Mrc1_Mac", "M14", "IFN_Mac", "Inf_Mac_Ly6cHi_Fn1Hi", "Inf_Mac_Ly6cLo_ACEHi", "cDC1", "cDC2", "pDC"))

	#Obj.markers$cluster <- factor(Obj.markers$cluster, levels = c("Trem2_hi_Mac", "Resident_Mac", "Inflammatory_Mac", "Mrc1_Mac", "M14", "IFN_Mac", "Inf_Mac_Ly6cHi_Fn1Hi", "Inf_Mac_Ly6cLo_ACEHi", "cDC1", "cDC2", "pDC"))

        top10 <- Obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

	#top10 <- top10[with(top10, order(cluster, -p_val)), ]

	#top10 <- read.delim(file = "cluster.lst.ann_top10Markers.xls.ann", row.names = 1, as.is = T)

        write.table(top10, file = paste0(AnnF, "_top10Markers.xls.ann"), col.names = NA, quote = FALSE, sep = "\t")

	#DefaultAssay(Obj) <- "integrated"

	all.genes <- rownames(Obj)

	Obj <- ScaleData(Obj, features = all.genes)

	#pp <- DoHeatmap(Obj, features = top10$gene, size = 2, raster = FALSE, group.colors = Colors[1:length(unique(Obj$category_anno_meta))]) + NoLegend() + theme(axis.text.y = element_text(size = 5)) + scale_fill_gradientn(colors = brewer.pal(n = 11, name = "RdBu")[c(9, 6, 3)])

        #pp <- DoHeatmap(Obj, features = top10$gene, size = 2, raster = FALSE, group.colors = Colors[1:length(unique(Obj$category_anno_meta))]) + NoLegend() + theme(axis.text.y = element_text(size = 5)) + scale_fill_gradientn(colors = brewer.pal(n = 11, name = "RdBu")[c(11:1)])

	pp <- DoHeatmap(Obj, features = top10$gene, size = 2, raster = FALSE, group.colors = Colors) + theme(axis.text.y = element_text(size = 5, face = "bold")) + scale_fill_gradientn(colors = brewer.pal(n = 11, name = "RdBu")[c(11:1)])

        #pp <- DoHeatmap(Obj, features = top10$gene, size = 2, raster = FALSE, group.colors = Colors[1:length(unique(Obj$category_anno_meta))]) + NoLegend() + theme(axis.text.y = element_text(size = 5)) + scale_fill_gradientn(colors = brewer.pal(n = 9, name = "YlOrRd")[c(1:9)])
	
	ggsave(paste0(AnnF, "_Markers_heatmap.ann.pdf"), pp, height = 9, width = 12)

        #DefaultAssay(Obj) <- "RNA"

}

#Expr <- GetAssayData(Obj, slot = "data")

Colors <- Colors[1:length(levels(Obj$category_anno_meta))]

names(Colors) <- levels(Obj$category_anno_meta)

pdf(paste0(AnnF, "_TSNEPlot_Ann.pdf"), width = 10)

	DimPlot(Obj, reduction = "umap", label = TRUE, repel = T, label.size = 7, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(Obj$category_anno_meta)]) + NoLegend() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()


pdf(paste0(AnnF, "_TSNEPlot_Blank.pdf"), width = 10)

        DimPlot(Obj, reduction = "umap", label = F, repel = T, label.size = 8, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(Obj$category_anno_meta)]) + NoLegend() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()

plist <- list()

for(ss in unique(Obj$orig.ident)){

	ObjTmp <- subset(Obj, cells = names(Obj$orig.ident[Obj$orig.ident == ss]))

	pp  <- DimPlot(ObjTmp, reduction = "umap", label = TRUE, repel = T, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(ObjTmp$category_anno_meta)]) + NoLegend() + ggtitle(ss)

	plist[[ss]] <- pp	

}

ggsave(file = paste0(AnnF, "_TSNEPlot_Sam.pdf"), arrangeGrob(grobs = plist, ncol = 2), width = 12, height = 4 * (as.integer(length(plist)/2)+1), limitsize = FALSE)


plist <- list()

for(ss in unique(Obj$SamInfo)){

        ObjTmp <- subset(Obj, cells = names(Obj$SamInfo[Obj$SamInfo== ss]))

        pp  <- DimPlot(ObjTmp, reduction = "umap", label = TRUE, repel = T, label.size = 6, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(ObjTmp$category_anno_meta)]) + NoLegend() + ggtitle(ss)

        plist[[ss]] <- pp

}

ggsave(file = paste0(AnnF, "_TSNEPlot_Condition.pdf"), arrangeGrob(grobs = plist, ncol = 2), width = 12, height = 4 * (as.integer(length(plist)/2)+1), limitsize = FALSE)


pdf(paste0(AnnF, "_TSNEPlot_Condition_mix.pdf"), width = 10)

        DimPlot(Obj, reduction = "umap", label = TRUE, repel = T, label.size = 4, pt.size = 0.3, group.by = "SamInfo", cols = Colors) + NoLegend() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))

dev.off()


plist <- list()

RamNum <- min(table(Obj$SamInfo))

for(ss in unique(Obj$SamInfo)){

        ObjTmp <- subset(Obj, cells = sample(names(Obj$SamInfo[Obj$SamInfo== ss]), RamNum))

        pp  <- DimPlot(ObjTmp, reduction = "umap", label = TRUE, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(ObjTmp$category_anno_meta)], repel = T) + NoLegend() + ggtitle(ss)

        plist[[ss]] <- pp

}

ggsave(file = paste0(AnnF, "_TSNEPlot_Condition.random.pdf"), arrangeGrob(grobs = plist, ncol = 5), width = 22, height = 3 * (as.integer(length(plist)/5)+1), limitsize = FALSE)


plist <- list()

RamNum <- min(table(Obj$SamInfo))

for(ss in unique(Obj$SamInfo)){

        ObjTmp <- subset(Obj, cells = sample(names(Obj$SamInfo[Obj$SamInfo== ss]), RamNum))

        pp  <- DimPlot(ObjTmp, reduction = "umap", label = F, pt.size = 0.3, group.by = "category_anno_meta", cols = Colors[levels(ObjTmp$category_anno_meta)], repel = T) + NoLegend() + ggtitle(ss)

        plist[[ss]] <- pp

}

ggsave(file = paste0(AnnF, "_TSNEPlot_Condition.random.blank.pdf"), arrangeGrob(grobs = plist, ncol = 5), width = 22, height = 3 * (as.integer(length(plist)/5)+1), limitsize = FALSE)


#ggsave(file = paste0(AnnF, "_TSNEPlot_Condition.random.pdf"), arrangeGrob(grobs = plist, ncol = 2), width = 12, height = 4 * (as.integer(length(plist)/2)))


FeaQC <- c("log10_nGene", "log10_nUMI", "percent.mt", "percent.rib")

        plist <- FeaturePlot(object = Obj, features = FeaQC, cols = c("grey", brewer.pal(n = 4, name = "Set3")[4]), reduction = "umap", pt.size = 0.3, combine = FALSE)

        ggsave(file = paste0(AnnF, "_Feature_plotQC.pdf"), arrangeGrob(grobs = plist, ncol = 3), width = 12, height = 2.3 * (as.integer(length(FeaQC)/3)+1), limitsize = FALSE)


        plist <- VlnPlot(object = Obj, features = FeaQC, pt.size = 0, combine = FALSE)

	plist <- lapply(plist, function(p) p <- p + theme(legend.position = "null"))

        ggsave(file = paste0(AnnF, "_Vln_plotQC.pdf"), arrangeGrob(grobs = plist, ncol = 1), width = 12, height = 3 * length(FeaQC), limitsize = FALSE)


	Idents(Obj) <- factor(paste(Obj$SamInfo, Obj$category_anno_meta, sep = "_"), levels = unlist(lapply(unique(Obj$category_anno_meta), function(x) paste(levels(Obj$SamInfo), x, sep = "_"))))

        plist <- VlnPlot(object = Obj, features = FeaQC, pt.size = 0, combine = FALSE)

        plist <- lapply(plist, function(p) p <- p + theme(legend.position = "null"))

        ggsave(file = paste0(AnnF, "_Vln_plotQC2.pdf"), arrangeGrob(grobs = plist, ncol = 1), width = length(levels(Idents(Obj))) * 0.5, height = 3 * length(FeaQC), limitsize = FALSE)

        #markers_c <- read.table(file = "top10Markers.xls", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

        #markers_cc <- intersect(markers_c$gene, rownames(Obj))

        #plist <- FeaturePlot(object = Obj, features = markers_cc, split.by = "orig.ident", cols = c("grey", brewer.pal(n = 4, name = "Set3")[4]), reduction = "umap", pt.size = 0.3, combine = FALSE)

        #ggsave(file = paste0(Args[3], "_Feature_plot.pdf"), arrangeGrob(grobs = plist, ncol = length(markers_cc)), height = length(unique(Obj$orig.ident)) * 4, width = 5 * length(markers_cc), limitsize = FALSE)

        #plist <- VlnPlot(object = Obj, features = markers_cc, pt.size = 0, combine = FALSE)

        #ggsave(file = paste0(Args[3], "_Vln_plot.pdf"), arrangeGrob(grobs = plist, ncol = 1), width = 12, height = 2.3 * length(markers_cc), limitsize = FALSE)


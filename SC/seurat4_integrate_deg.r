
#! R/4.0.3

library(Seurat)
library(RColorBrewer)
library(lattice)
library(ggplot2)
library(gridExtra)
library(cowplot)

Args <- commandArgs(trailingOnly = TRUE)

options(future.globals.maxSize = 1000 * 1024^2)

ObjF <- Args[1]
SamGroupF <- Args[2]
CompareS <- Args[3]
CellsS <- Args[4]

Obj <- readRDS(file = ObjF)

dim(Obj)

SamInfo <- read.delim(file = SamGroupF, row.names = 1, stringsAsFactors = FALSE)

Obj <- subset(Obj, subset = orig.ident %in% rownames(SamInfo))

dim(Obj)

Obj$SamInfo <- SamInfo[Obj$orig.ident, ]

Compare <- unlist(strsplit(CompareS, split = ";"))

#Obj <- subset(Obj, subset = SamInfo == Compare)

Obj <- subset(Obj, subset = SamInfo %in% Compare)

table(Obj$SamInfo)

Obj$celltype.CC <- paste(Idents(Obj), Obj$SamInfo, sep = "_")

Obj$celltype <- Idents(Obj)

Idents(Obj) <- "celltype.CC"

CountId <- table(Idents(Obj))

DefaultAssay(Obj) <- "RNA"

DEG <- list()

#CCs <- unlist(strsplit(Args[2], split = ";"))

Cells <- unlist(strsplit(CellsS, split = ";"))

for (cc in Cells){

        print(cc)

	if(!is.na(CountId[paste(cc, Compare[1], sep = "_")]) & !is.na(CountId[paste(cc, Compare[2], sep = "_")])){

	#if(CountId[paste(cc, Compare[1], sep = "_")] > 3 & CountId[paste(cc, Compare[2], sep = "_")] > 3){
        if(CountId[paste(cc, Compare[1], sep = "_")] > 20 & CountId[paste(cc, Compare[2], sep = "_")] > 20){

        	DEG[[cc]] <- FindMarkers(Obj, ident.1 = paste(cc, Compare[1], sep = "_"), ident.2 = paste(cc, Compare[2], sep = "_"), verbose = FALSE, logfc.threshold = 0, min.pct = 0.1, slot = "data")

	}

	}

}

saveRDS(DEG, file = paste0(Compare[1], "_", Compare[2], "_DEG.rds"))

for(cc in names(DEG)){

        write.table(DEG[[cc]], file = paste0(Compare[1], "_", Compare[2], "_", cc, "_DEG.xls"), sep = "\t", quote = FALSE, col.names = NA)

}

DEGsig <- lapply(DEG, function(x) rownames(x[x$avg_logFC >= 0.25, ])[1:5])

DEGsig <- lapply(DEGsig, function(x) x[!(is.na(x))])


Idents(Obj) <- Obj$celltype

pdf(paste0(Compare[1], "_", Compare[2], "_Dot_plot_DEG.pdf"), width = 20)

        DotPlot(Obj, features = unique(unlist(DEGsig)), cols = brewer.pal(12, "Set3")[c(4, 5)], dot.scale = 8, split.by = "SamInfo") + RotatedAxis()

dev.off()

Plots <- list()

for(cc in names(DEGsig)){

        if(length(DEGsig[[cc]])==0){DEGsig[[cc]] <- NULL}

}

for(cc in names(DEGsig)){

        theme_set(theme_cowplot())

        CD3Cell <- subset(Obj, idents = cc)

        Idents(CD3Cell) <- "SamInfo"

        avg.CD3Cell <- log1p(AverageExpression(CD3Cell, verbose = FALSE)$RNA)

        avg.CD3Cell$gene <- rownames(avg.CD3Cell)

        genes.to.label = DEGsig[[cc]]

        p1 <- ggplot(avg.CD3Cell, aes_string(colnames(avg.CD3Cell)[1], colnames(avg.CD3Cell)[2])) + geom_point(color = brewer.pal(12, "Set3")[4]) + ggtitle(cc)

        p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, size = 5)

        Plots[[cc]] <- p1

        #ggsave(paste0("Diff_plot_", cc, ".pdf"), p1)

}

pdf(paste0(Compare[1], "_", Compare[2], "_Diff_plot_DEG.pdf"), width = 16, height = 3 * (as.integer(length(DEGsig)/4) + 1))

        CombinePlots(plots = Plots, ncol = 4)

dev.off()

pdf(paste0(Compare[1], "_", Compare[2], "_Vln_plot_DEG.pdf"), width = 12, height = 3 * length(unique(unlist(DEGsig))))

        plots <- VlnPlot(Obj, features = unique(unlist(DEGsig)), split.by = "SamInfo", group.by = "category_anno_meta", pt.size = 0, combine = FALSE)

        CombinePlots(plots = plots, ncol = 1)

dev.off()


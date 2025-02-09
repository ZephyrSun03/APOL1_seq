#! R/4.0.3

Args <- commandArgs(trailingOnly = TRUE)

FileList <- Args[1]
FileType <- Args[2]
#MetaData <- Args[3]
MarkerData <- Args[3]
Resolution <- Args[4] ## 0.8
DimNum <- Args[5]
DimNum2 <- Args[6]
MinCell <- Args[7]      ## 3
MinFeature <- Args[8]   ## 200
MaxFeature <- Args[9]   ## 2500
MTPerc <- Args[10] ## 15
Genomes <- Args[11]
VarFeatureN <- as.integer(Args[12]) ## 2000

library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(lattice)
library(gridExtra)
#library(harmony)

File.list <- read.table(file = FileList, header = FALSE, stringsAsFactors = FALSE, sep = "\t", row.names = 2)

QCPlots <- list()

QCPlots2 <- list()

Data.list <- list()


for(x in rownames(File.list)){

	print(x)

	if(File.list[x, 2] == "h5"){

		Expr <- Read10X_h5(file = File.list[x, 1])

	}else if(File.list[x, 2] == "RDS"){

		Expr <- readRDS(file = File.list[x, 1])

	}else if(File.list[x, 2] == "table"){

		Expr <- read.delim(file = File.list[x, 1], header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

	}else if(File.list[x, 2] == "csv"){

        	Expr <- read.csv(file = File.list[x, 1], header = TRUE, row.names = 1, check.names = FALSE)

	}else if(File.list[x, 2] == "10x"){

		Expr <- Read10X(data.dir = File.list[x, 1])

                if(any(grepl("Gene Expression", names(Expr)))){

                        Expr <- Expr[["Gene Expression"]]

                }


	}else if(File.list[x, 2] == "10xMulti"){

		inputdata.10x <- Read10X(data.dir = File.list[x, 1])

		Expr <- inputdata.10x$`Gene Expression`

	}

	if(Genomes != "All"){

		Genomes <- unlist(strsplit(Genomes, split = ";"))

		Expr <- Reduce(rbind, lapply(Genomes, function(x) {Tmp <- Expr[grep(x, rownames(Expr)), ]; rownames(Tmp) <- gsub(x, "", rownames(Tmp)); Tmp}))

	}

	colnames(Expr) <- paste(x, colnames(Expr), sep = "_")

	print("Raw Expr")

	print(dim(Expr))

	#Object <- CreateSeuratObject(counts = Expr, project = x, min.cells = as.integer(MinCell), min.features = as.integer(MinFeature), names.delim = ";")

	Object <- CreateSeuratObject(counts = Expr, project = x, min.cells = as.integer(MinCell), names.delim = ";")

	#Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^MT-") + PercentageFeatureSet(Object, pattern = "^mt-")

	Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^MT-") + PercentageFeatureSet(Object, pattern = "^mt-") + PercentageFeatureSet(Object, pattern = "^hg38--MT-")

	Object[["percent.rib"]] <- PercentageFeatureSet(Object, pattern = "^RPS") + PercentageFeatureSet(Object, pattern = "^RPL")+ PercentageFeatureSet(Object, pattern = "^Rps-") + PercentageFeatureSet(Object, pattern = "^Rpl") + PercentageFeatureSet(Object, pattern = "^hg38--RPS") + PercentageFeatureSet(Object, pattern = "^hg38--RPL")

	Object$log10_nGene <- log10(Object$nFeature_RNA + 1)

	Object$log10_nUMI <- log10(Object$nCount_RNA + 1)

	QCPlots <- c(QCPlots, VlnPlot(Object, features = c("log10_nGene", "log10_nUMI", "percent.mt", "percent.rib"), pt.size = 0, ncol = 4, combine = FALSE))
	
	Object <- subset(Object, subset = nFeature_RNA > as.integer(MinFeature) & nFeature_RNA < as.integer(MaxFeature) & percent.mt < as.numeric(MTPerc))

	print(paste(dim(Object)[2]," cells remained!", sep=""))

	print(paste(dim(Object)[1]," genes remained!", sep=""))

	QCPlots2 <- c(QCPlots2, VlnPlot(Object, features = c("log10_nGene", "log10_nUMI", "percent.mt", "percent.rib"), pt.size = 0, ncol = 4, combine = FALSE))

	#Object <- NormalizeData(Object, normalization.method = "LogNormalize", scale.factor = 10000)
	
        #Object <- FindVariableFeatures(Object, selection.method = "vst", nfeatures = 2000)

	Data.list[[x]] <- Expr

}

ggsave(file = "1.QC_plot.pdf", arrangeGrob(grobs = QCPlots, ncol = 4), width = 4*4, height = 5 * nrow(File.list), limitsize = FALSE)

ggsave(file = "1.QC_plot2.pdf", arrangeGrob(grobs = QCPlots2, ncol = 4), width = 4*4, height = 5 * nrow(File.list), limitsize = FALSE)

names(Data.list) <- rownames(File.list)

Group <- c()

for(i in rownames(File.list)){

	#Anno <- unlist(read.delim(file = File.list[i, "Ann"], row.names = 1, stringsAsFactors = FALSE))

	#names(Anno) <- levels(Data.list[[i]])

	#pbmc <- RenameIdents(Data.list[[i]], Anno)

	tmp <- rep(i, ncol(Data.list[[i]]))

	names(tmp) <- colnames(Data.list[[i]])

	Group <- c(Group, tmp)

}

#for (i in 1:length(Data.list)) {
#
#	Data.list[[i]] <- NormalizeData(Data.list[[i]], verbose = FALSE)
#
#	Data.list[[i]] <- FindVariableFeatures(Data.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#
#}

saveRDS(Data.list, file = "Data.list.rds")

#Genes <- unique(unlist(lapply(Data.list, function(x) rownames(x))))

Genes <- Reduce(intersect, lapply(Data.list, rownames))

Data.list <- lapply(Data.list, function(x) x[Genes, ])

ExprAll <- Reduce(cbind, Data.list)

#names(Group) <- colnames(ExprAll)

print(paste0("Original:", nrow(ExprAll), " Genes ", ncol(ExprAll), "Cells!"))

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes



Object <- CreateSeuratObject(counts = ExprAll, project = "PBMC", min.cells = as.integer(MinCell), names.delim = ";") ## updated in V3

names(Group) <- colnames(Object)

Object$orig.ident <- Group

str(Object)

        Object$log10_nGene <- log10(Object$nFeature_RNA + 1)

        Object$log10_nUMI <- log10(Object$nCount_RNA + 1)

Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = "^MT-") + PercentageFeatureSet(Object, pattern = "^mt-") + PercentageFeatureSet(Object, pattern = "^hg38--MT-")

Object[["percent.rib"]] <- PercentageFeatureSet(Object, pattern = "^RPS") + PercentageFeatureSet(Object, pattern = "^RPL")+ PercentageFeatureSet(Object, pattern = "^Rps") + PercentageFeatureSet(Object, pattern = "^Rpl") + PercentageFeatureSet(Object, pattern = "^hg38--RPS") + PercentageFeatureSet(Object, pattern = "^hg38--RPL")

pdf("1.QC.pdf", width = 16)

        #VlnPlot(object = Object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, x.lab.rot = TRUE)

        VlnPlot(Object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"), ncol = 4, pt.size = 0)

dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

pdf("2.QC_Cor.pdf")

        plot1 <- FeatureScatter(Object, feature1 = "nCount_RNA", feature2 = "percent.mt")
        plot2 <- FeatureScatter(Object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        CombinePlots(plots = list(plot1, plot2))

dev.off()


Object <- subset(Object, subset = nFeature_RNA > as.integer(MinFeature) & nFeature_RNA < as.integer(MaxFeature) & percent.mt < as.numeric(MTPerc))

print(paste0("QCed:", nrow(Object), " Genes ", ncol(Object), "Cells!"))

dim(Object)

print("QC done!")

Object <- NormalizeData(Object, normalization.method = "LogNormalize", scale.factor = 10000)

print("Normalization done!")


pdf("3.variant_gene.pdf", width = 14)

        #Object <- FindVariableGenes(object = Object, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 6, y.cutoff = 1)

        Object <- FindVariableFeatures(Object, selection.method = "vst", nfeatures = VarFeatureN)

        #if(VarialMan != "No"){

        #        Vars <- read.table(file = VarialMan, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

        #        print(paste0(length(intersect(Vars$Gene, rownames(Object))), "Vars found!"))

        #        VariableFeatures(Object) <- intersect(Vars$Gene, rownames(Object))

        #}

        # Identify the 10 most highly variable genes

        top10 <- head(VariableFeatures(Object), 10)

        # plot variable features with and without labels

        plot1 <- VariableFeaturePlot(Object)

        LabelPoints(plot = plot1, points = top10, repel = TRUE)

        #CombinePlots(plots = list(plot1, plot2))

dev.off()


all.genes <- rownames(Object)

Object <- ScaleData(Object, features = all.genes)

print("Scaling done!")


Object <- RunPCA(Object, features = VariableFeatures(object = Object))

# Examine and visualize PCA results a few different ways
print(Object[["pca"]], dims = 1:5, nfeatures = 5)


pdf("4.VizPCA.pdf", width = 12, height = 12)

        #VizPCA(object = Object, pcs.use = 1:2)

        VizDimLoadings(Object, dims = 1:10, reduction = "pca")

dev.off()

pdf("5.PCAPlot.pdf")

        #PCAPlot(object = Object, dim.1 = 1, dim.2 = 2)

        DimPlot(Object, reduction = "pca")

dev.off()

#pdf("5.Harmony_Coverage.pdf")
#
#Object <- Object %>% RunHarmony("orig.ident", plot_convergence = TRUE)
#
#dev.off()

#harmony_embeddings <- Embeddings(Object, 'harmony')

#harmony_embeddings[1:5, 1:5]

#pdf("5.Harmony.pdf")
#
#	p1 <- DimPlot(object = Object, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
#	p2 <- VlnPlot(object = Object, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
#
#plot_grid(p1,p2)
#
#dev.off()


pdf("6.PCHeatmap.pdf", height = 14)

        #PCHeatmap(object = Object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

        DimHeatmap(Object, reduction = "pca", dims = 1:20, cells = 500, balanced = TRUE)

dev.off()


#Object <- JackStraw(Object, num.replicate = 100)
#Object <- ScoreJackStraw(Object, dims = 1:20)
#
#pdf("7.JackStraw.pdf")
#
#        JackStrawPlot(Object, dims = 1:20)
#
#dev.off()

pdf("8.PCElbowPlot.pdf")

        #PCElbowPlot(object = Object)

        ElbowPlot(Object, reduction = "pca")

dev.off()

#print("Harmony done!")

Object <- FindNeighbors(Object, reduction = "pca", dims = 1:as.integer(DimNum2))

Object <- FindClusters(Object, resolution = as.numeric(Resolution))

print(paste(nlevels(Idents(Object)), "clusters was identified!"))

write.table(Idents(Object), file = "cluster.lst", sep = "\t", quote = F, col.names = F)

print("Clustering done!")


Object <- RunUMAP(Object, reduction = "pca", dims = 1:as.numeric(DimNum2))

#save(Object, file = "Object.rda")

pdf(paste0("9.TSNEPlot_resolution", Resolution, "_Dim", DimNum2, ".pdf"), width = 12)

       DimPlot(Object, reduction = "umap", label = TRUE)

dev.off()

Colors <- c(brewer.pal(8, "Set2"), brewer.pal(9, "Set1"), brewer.pal(12, "Set3"))

pdf(paste0("10.TSNEPlot_resolution_ori", Resolution, "_Dim", DimNum2, ".pdf"), width = 12)

	DimPlot(Object, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = FALSE, cols = Colors[1:(length(unique(Object$orig.ident)))])

dev.off()


        #DefaultAssay(Object) <- "RNA"

        markers_c <- read.table(file = MarkerData, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

        markers_cc <- intersect(markers_c$Gene, rownames(Object))

	plist <- FeaturePlot(object = Object, features = markers_cc, cols = c("grey", brewer.pal(n = 4, name = "Set3")[4]), reduction = "umap", pt.size = 0.3, combine = FALSE)
	
	ggsave(file = "11.Feature_plot.pdf", arrangeGrob(grobs = plist, ncol = 3), width = 12, height = 2.3 * (as.integer(length(markers_cc)/3)+1), limitsize = FALSE)

        plist <- VlnPlot(object = Object, features = markers_cc, pt.size = 0, combine = FALSE)

        ggsave(file = "12.Vln_plot.pdf", arrangeGrob(grobs = plist, ncol = 1), width = 12, height = 2.3 * length(markers_cc), limitsize = FALSE)
	
        #DefaultAssay(Object) <- "integrated"

#        Object.markers <- FindAllMarkers(Object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#
#        write.table(Object.markers, file = "Markers.lst", sep = "\t", quote = F, col.names = NA)
#
#        print(paste(dim(Object.markers)[1], "markers found!"));
#
#        top10 <- Object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#
#        write.table(top10, file = "top10Markers.xls", col.names = NA, quote = FALSE, sep = "\t")

saveRDS(Object, file = "Obj.rds")

#pdf("13.Markers_heatmap.pdf")
#
#        DoHeatmap(Object, features = top10$gene, size = 2, raster = FALSE) + NoLegend() + scale_fill_gradientn(colors = brewer.pal(n = 11, name = "RdBu")[c(9, 6, 3)]) + theme(axis.text.y = element_text(size = 3))
#
#dev.off()

#save(Object, file = "Object.rda")

#saveRDS(Object, file = "integrated.rds")



#!/usr/bin/R

library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(ggrepel)
#library(gg3D)
#library(plotly)
#library(DESeq2)

Args <- commandArgs(trailingOnly = TRUE)

ExprF <- Args[1]
ThresCPM <- as.numeric(Args[2]) ## 1
CheckS <- Args[3]
GenesF <- Args[4]
FactorF <- Args[5]
KeyS <- Args[6]

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", header = F, row.names = 1)

Expr <- read.table(file = ExprF, sep = "\t", header = TRUE, row.names = 1, check.names = F)

dim(Expr)

dge <- DGEList(counts = Expr)

dge <- calcNormFactors(dge)

Expr.norm <- cpm(dge, log = T, normalized.lib.sizes = T)

drop <- which(apply(cpm(dge), 1, max) >= ThresCPM)

length(drop)

dge <- dge[drop,]

if(GenesF != "All"){

	GeneLis <- unlist(read.table(file = GenesF, as.is = T))

	length(intersect(GeneLis, rownames(dge)))

	dge <- dge[rownames(dge) %in% GeneLis, ]
}

#keep <- rowSums(counts(dge)) >= ThresReadCount

#dge <- dge[keep,]


Expr.norm <- Expr.norm[rownames(dge), ] 

#write.table(Expr.norm, file = "Expr.cpm.xls", sep = "\t", quote = F, col.names = NA)

## PCA analysis

Factor <- read.table(file = FactorF, sep = "\t", header = TRUE, row.names = 1, as.is = T)

#row.names(Factor) <- Factor$Sam

int_PCA <- prcomp(t(Expr.norm))

int_PCA_x <- merge(Factor, int_PCA$x[, 1:4], by.x = 0, by.y = 0)

Variance <- summary(int_PCA)$importance["Proportion of Variance", ] * 100

head(Variance)

#checks <- c("Cell", "Stiffness")

checks <- unlist(strsplit(CheckS, split = ";"))

for(check in checks){

        p<-ggplot(int_PCA_x,aes_string(x="PC1",y="PC2",color=check, label="Row.names"))

        p<-p+geom_point(size = 5) + geom_text_repel(size=3) + theme_minimal() + xlab(paste0("PC1 ", Variance[1], "%")) + ylab(paste0("PC2 ", Variance[2], "%"))

        ggsave(paste0(check, "_", KeyS, "_PCA_dim_1_2.pdf"), p, width = 4, height = 4)

}

#for(check in checks){
#
#	fig <- plot_ly(int_PCA_x, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Cell, colors = Colors[c("Red", "Blue", "Yellow"), ])
#
#	#fig <- fig %>% add_markers()
#
#	fig <- fig %>% layout(scene = list(xaxis = list(title = paste0("PC1 ", Variance[1], "%")),
#                     yaxis = list(title = paste0("PC2 ", Variance[2], "%")),
#                     zaxis = list(title = paste0("PC3 ", Variance[3], "%"))))
#
#	orca(fig, paste0(check, "_PCA_dim_1_2_3.pdf"))
#
#}

## MDS

mds <- plotMDS(dge)

        for(nn in checks){

                Tab <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = Factor[names(mds$x), nn])

                pp <- ggplot(Tab, aes(x = Dim1, y = Dim2, color = Group)) +

                        geom_point() +

			xlab("MDS1") + ylab("MDS2") +

                        geom_text_repel(aes(label = rownames(Tab)), size = 3) +

                        theme_minimal()

                ggsave(paste0(KeyS, "_", nn, "_MDSPlot.pdf"))

        }


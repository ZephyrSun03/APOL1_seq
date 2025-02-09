
Args <- commandArgs(trailingOnly = T)

library("limma")
library("variancePartition")
library(edgeR)
library(ggrepel)

ExprF <- Args[1]
FacF <- Args[2]
ThresCPMN <- as.numeric(Args[3])
FocusS <- Args[4]
CheckS <- Args[5]

Colors <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/ColorLibrary", row.names = 1, header = F, as.is = T)

Expr <- read.delim(file = ExprF, row.names = 1, as.is = T, check.names = F)
Fac <- read.delim(file = FacF, row.names = 1, as.is = T)

Expr <- Expr[, rownames(Fac)]

dge <- DGEList(counts = Expr)

dge <- calcNormFactors(dge)

drop <- which(apply(cpm(dge), 1, max) >= ThresCPMN)

dge <- dge[drop,]

print("Filtered")

dim(dge) # number of genes left

design <- model.matrix(~ 0 + Fac[[FocusS]])

pdf("Voom.pdf")

	dge2 <- voom(dge, design, plot = T)

dev.off()


form <- ~ (1 | Condition) + (1 | Batch)

varPart <- fitExtractVarPartModel(dge2$E, form, Fac)


# sort variables (i.e. columns) by median fraction
#       of variance explained
vp <- sortCols(varPart)

# Figure 1a
# Bar plot of variance fractions for the first 10 genes

pp <- plotPercentBars(vp[1:10, ])

ggsave("GeneVar.pdf", pp)


# Figure 1b
# violin plot of contribution of each variable to total variance

pp <- plotVarPart(vp)

ggsave("VariableVar.pdf", pp)


# subtract out effect of Batch

fit <- lmFit(dge2, model.matrix(~Batch, Fac))

res <- residuals(fit, dge2)

write.table(res, file = paste0(ExprF, "_Brm.xls"), sep = "\t", quote = F, col.names = NA)

# fit model on residuals

form <- ~ (1 | Condition) + (1 | Batch)

varPartResid <- fitExtractVarPartModel(res, form, Fac)

vp <- sortCols(varPartResid)

pp <- plotPercentBars(vp[1:10, ])

ggsave("GeneVar.Brm.pdf", pp)

pp <- plotVarPart(vp)

ggsave("VariableVar.Brm.pdf", pp)


int_PCA <- prcomp(t(res))

int_PCA_x <- merge(Fac, int_PCA$x[, 1:4], by.x = 0, by.y = 0)

Variance <- summary(int_PCA)$importance["Proportion of Variance", ] * 100

head(Variance)

#checks <- c("Cell", "Stiffness")

checks <- unlist(strsplit(CheckS, split = ";"))

for(check in checks){

        p<-ggplot(int_PCA_x,aes_string(x="PC1",y="PC2",color=check, label="Row.names"))

        p<-p + geom_point(size = 5) + geom_text_repel(size=3) + theme_minimal() + xlab(paste0("PC1 ", Variance[1], "%")) + ylab(paste0("PC2 ", Variance[2], "%")) + theme_classic() + scale_color_manual(values = Colors[c("Blue", "Red", "SteelDark"), ])

        ggsave(paste0(check, "_PCA_dim_1_2.Brm.pdf"), p, width = 4, height = 4)

}


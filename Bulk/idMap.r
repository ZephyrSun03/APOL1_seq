
Expr <- read.table(file = "Expr.xls", row.names = 1, header = TRUE, sep = "\t", check.names = F)

quantile(colSums(Expr))

Ann <- read.table(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/gene_symbol.txt", header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, row.names = 1)

#Ann <- read.delim(file = "/sc/arion/projects/zhangw09a/PANDA/db_ZS/Ensembl/mm39/Mus_musculus.GRCm39.105.gff3.corr.gff3.map.ann", as.is = T, row.names = 2)

#Expr <- merge(Ann, Expr, by.x =1 , by.y = 0)

Symbol <- Ann[rownames(Expr), "HGNC.symbol"]

length(which(Symbol == ""))

length(which(is.na(Symbol)))

length(which(duplicated(Symbol)))


Expr <- Expr[which(Symbol != ""), ]

Symbol <- Symbol[which(Symbol != "")]

Expr <- Expr[!duplicated(Symbol), ]

Symbol <- Symbol[!duplicated(Symbol)]

Expr <- Expr[!is.na(Symbol), ]

Symbol <- Symbol[!is.na(Symbol)]

rownames(Expr) <- Symbol

#row.names(Expr) <- make.names(Expr[[2]], unique = TRUE)

#Expr <- Expr[, c(-1, -2)]

write.table(Expr, file = "Expr.xls.sym", sep = "\t", col.names = NA, quote = FALSE)

library(edgeR)

dge <- DGEList(counts = Expr)

#dge <- calcNormFactors(dge, method = "none")

Expr.norm <- cpm(dge, log = T, normalized.lib.sizes = T)

write.table(Expr.norm, file = "Expr.xls.sym.logcpm.xls", sep = "\t", col.names = NA, quote = FALSE)


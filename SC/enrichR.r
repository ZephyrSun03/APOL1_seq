
# R/4.0.3

Args <- commandArgs(trailingOnly = T)

library(enrichR)

#DBs <- listEnrichrDbs()

#Genes <- read.delim(file = Args[1], row.names = 1, header = F, as.is = T)

Genes <- read.delim(file = Args[1], header = F, as.is = T)

Genes <- Genes[!is.na(Genes[[1]]), ]

Genes <- Genes[!duplicated(Genes[[1]]), ]

rownames(Genes) <- Genes[[1]]

Genes <- Genes[, -1, drop = F]

#dbs <- c("GO_Biological_Process_2021", "WikiPathway_2021_Human", "KEGG_2021_Human", "MSigDB_Hallmark_2020", "BioCarta_2016", "Reactome_2016", "")

#dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "WikiPathway_2021_Human", "MSigDB_Hallmark_2020", "KEGG_2021_Human", "BioCarta_2016", "Reactome_2016", "COVID-19_Related_Gene_Sets_2021", "COVID-19_Related_Gene_Sets", "Cancer_Cell_Line_Encyclopedia", "NCI-60_Cancer_Cell_Lines", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "NCI-Nature_2016", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019")

dbs <- c("GO_Biological_Process_2021", "GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "WikiPathway_2021_Human", "MSigDB_Hallmark_2020", "KEGG_2021_Human", "BioCarta_2016", "Reactome_2016", "BioPlanet_2019", "Cancer_Cell_Line_Encyclopedia", "NCI-60_Cancer_Cell_Lines", "GTEx_Tissue_Sample_Gene_Expression_Profiles_down", "GTEx_Tissue_Sample_Gene_Expression_Profiles_up", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "NCI-Nature_2016", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019")

#dbs <- c("WikiPathway_2021_Human", "WikiPathways_2019_Human", "MSigDB_Hallmark_2020", "WikiPathways_2016", "KEGG_2021_Human", "KEGG_2019_Human", "KEGG_2016", "KEGG_2013", "BioCarta_2016", "BioCarta_2015", "Reactome_2016", "Reactome_2013", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019")

enriched <- enrichr(rownames(Genes), dbs)

#saveRDS(enriched, file = "enriched.rds")

Temp <- unlist(lapply(enriched, function(x) nrow(x[x$P.value <= 0.1, ])))

for(nn in names(Temp)[Temp == 0]){

	enriched[[nn]] <- NULL

}

#Sigs <- lapply(names(enriched), function(x) data.frame(enriched[[x]][enriched[[x]][["P.value"]] <= 0.1, ], DB = x))

Sigs <- Reduce(rbind, lapply(names(enriched), function(x) data.frame(enriched[[x]][enriched[[x]][["P.value"]] <= 0.1, ], DB = x)))

DEG <- Genes

Sigs$Pos <- sapply(Sigs$Genes, function(x) {Genes <- unlist(strsplit(x, split = ";")); length(which(DEG[Genes, ]>0))})

Sigs$Neg <- sapply(Sigs$Genes, function(x) {Genes <- unlist(strsplit(x, split = ";")); length(which(DEG[Genes, ]<0))})

Sigs$Count <- Sigs$Pos + Sigs$Neg

write.table(Sigs, file = paste0(Args[1], "_EnrichR.sig.xls"), sep = "\t", quote = F, row.names = F)


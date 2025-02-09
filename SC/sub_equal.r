
## R/4.0.3

#Args <- commandArgs(trailingOnly = T)

Args <- c("/sc/arion/scratch/sunz04/APOL1/SCNew/BootStrap/cluster.lst.ann_Obj_post.rds.hto.rds", "/sc/arion/projects/GOCAR/Sun/5.DRWork/6.APOL1/SCNew/1.Cluster/Sample.info.HTO")

library(Seurat)

ObjF <- Args[1]
SamGroupF <- Args[2]
#KeyS <- Args[3]

#set.seed(as.integer(KeyS))

Obj <- readRDS(file = ObjF)

dim(Obj)

SamInfo <- read.delim(file = SamGroupF, row.names = 1, stringsAsFactors = FALSE)

Obj <- subset(Obj, subset = orig.ident %in% rownames(SamInfo))

dim(Obj)

Obj$SamInfo <- SamInfo[Obj$orig.ident, ]

Cells <- levels(Obj$category_anno_meta)

RList <- c()

for(nn in 1:100){

	#set.seed(nn)

	Randoms <- c()
	
	for(cc in Cells){
	
		Tmp <- Obj@meta.data[Obj$category_anno_meta == cc, ]
	
		Counts <- table(Tmp$SamInfo)
	
		Randoms <- c(Randoms, rownames(Tmp)[Tmp$SamInfo == names(Counts)[which(Counts == min(Counts))]], sample(rownames(Tmp)[Tmp$SamInfo == names(Counts)[which(Counts == max(Counts))]], min(Counts)))
	
	}

	RList[[nn]] <- Randoms
	
	#Obj.sub <- subset(Obj, cells = Randoms)
	
	#table(Obj.sub@meta.data[, c("category_anno_meta", "SamInfo")])
	
	#saveRDS(Obj.sub, file = paste0(ObjF, "_Random", nn, ".rds"))

}

for(nn in 1:100){

	Obj.sub <- subset(Obj, cells = RList[[nn]])

        table(Obj.sub@meta.data[, c("category_anno_meta", "SamInfo")])

        saveRDS(Obj.sub, file = paste0(ObjF, "_Random", nn, ".rds"))
	
}

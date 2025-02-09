
Args <- commandArgs(trailingOnly = T)

ListF <- Args[1]

library(SoupX)
#library(DropletUtils)

Lists <- read.delim(file = ListF, row.names = 2, as.is = T, header = F)

for(nn in rownames(Lists)){

	tryCatch(

	{sc = load10X(Lists[nn, 1])
	
	sc = autoEstCont(sc)
	
	out = adjustCounts(sc)
	
	cntSoggy <- apply(sc$toc, 1, function(x) sum(x[x>0]))
	
	cntStrained <- apply(out, 1, function(x) sum(x[x>0]))
	
	Tab <- data.frame(Ori = cntSoggy, Cleaned = cntStrained)
	
	Tab$Rate <- (Tab$Ori - Tab$Cleaned)/Tab$Ori
	
	Tab <- Tab[order(Tab$Rate, decreasing = T), ]
	
	write.table(Tab, file = paste0(nn, "_stat.xls"), sep = "\t", quote = F, col.names = NA)
	
	saveRDS(out, file = paste0(nn, "_cleaned.rds"))

	}, error = function(e){cat("ERROR :",conditionMessage(e), "\n")}

	)

}

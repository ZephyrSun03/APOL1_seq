
Args <- commandArgs(trailingOnly = T)

File <- Args[1]

Lsts <- read.table(file = File, row.names = 2, as.is = T)

Stats <- lapply(rownames(Lsts), function(x) read.delim(file = paste0(x, "_stat.xls"), row.names = 1))

names(Stats) <- rownames(Lsts)

StatsTabs <- lapply(Stats, function(x) apply(x, 2, function(y) quantile(y, na.rm = T, probs = seq(0, 1, by = 0.01))))

for(nn in names(StatsTabs)){

	colnames(StatsTabs[[nn]]) <- paste0(colnames(StatsTabs[[nn]]), "_", nn)

}

StatsTab <- Reduce(cbind, StatsTabs)

write.table(StatsTab, file = "Ambient_stat.xls", sep = "\t", quote = F, col.names = NA)


lapply(Stats, function(x){x[head(order(x[["Rate"]], decreasing = T), 20), ]})


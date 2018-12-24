# Rscript make_output_table_betase.r 

args<-commandArgs(TRUE)
outpath <- args[1]
pval <- args[2]
prefix1 <- args[3]
prefix2 <- args[4]

errot <- "betase.norm" 

traits <- c("height")
pstatusj <- "_pruned"
pvals <- c(pval)
prefixes <- c(prefix1, prefix2)


for (traiti in traits){
	for (pvalk in pvals){
		
		fname <- paste(outpath, prefixes[1], ".",errot, ".pval.", pvalk, sep="")
		df <- read.table(fname, header=TRUE)
		outdf <-  data.frame(df$populations)
		
		
		for (p in prefixes){
			fname <- paste(outpath, p, "." , errot, ".pval.", pvalk, sep="")
			df <- read.table(fname, header=TRUE)
			outdf <- cbind(outdf, df$scores, df$lower, df$upper)
			write.table(outdf, file=paste(outpath, traiti, pstatusj, ".", errot,".summary.pval.", pvalk, sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
		}
	}
}


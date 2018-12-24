# Rscript visualize_prs_betase.r 


### assign file and path names ###
args<-commandArgs(TRUE)
mpath <- args[1]
outpath <- args[2]
idsf <- args[3]
prefix <- args[4]
pval <- args[5]
outtype <- args[6]
x <- c("CEU", "GBR", "IBS", "TSI")


### read in individual assignment file ###
ids <- read.table(idsf, header=TRUE)	
npops <- length(unique(ids$assignment))


### obtain original scores 
fname2 <- paste(prefix,".effectacounts.indsubset.subpops.gwas.pvalthres.", pval,".prs", sep="")
f2 <- paste(mpath, fname2, sep="")

### Subpopulation analysis ###
## read in file with subpopulation counts 
fields <- count.fields(f2, sep = '\t')
maxf <- max(fields)
df2 <- read.table(f2, header=FALSE, col.names = paste0("V",seq_len(maxf)), fill = TRUE)

## pull number and names of subpopulations
n <- df2[2,]
labels <- n[n != ""]
nlabels <- length(labels)
sumso <- as.numeric(t(df2[3,1:nlabels]))
totso <- as.numeric(t(df2[3,(nlabels+1):(nlabels*2)]))
countso <- as.numeric(t(df2[3,((nlabels*2)+1):(nlabels*3)]))
meano <- as.numeric(t(df2[3,((nlabels*3)+1):((nlabels*3)+2)]))[1]
Vao <- as.numeric(t(df2[3,((nlabels*3)+1):((nlabels*3)+2)]))[2]
varZ <- as.numeric(t(df2[3,((nlabels*3)+3):((nlabels*3)+3+nlabels-1)]))
normprso <- (sumso-meano)/sqrt(Vao)
se <- sqrt(varZ) 
normcilow <- ((sumso-1.96*se) - meano)/sqrt(Vao)
normcihi <- ((sumso + 1.96*se) - meano)/sqrt(Vao)
 

## compute times for each subpopulation
avg_dates <- c()
counter <- 1
for (labeli in labels){
	dates = ids$date[ids$assignment==labeli]
	avg_dates[counter] <- mean(dates)
	counter <- counter + 1
}
avg_dates <- avg_dates * -1

## create data frame for subpopulations

inds <- match(x, labels)


outdf <- data.frame(labels[inds], sumso[inds], se[inds], sumso[inds]-1.96*se[inds], sumso[inds]+1.96*se[inds],avg_dates[inds])
names(outdf) <- c("populations", "scores", "se", "lower", "upper","times")
write.table(outdf, file=paste(outpath, outtype, ".betase.pval.", pval, sep=""), quote=FALSE, row.names=FALSE)

outdf2 <- data.frame(labels[inds], normprso[inds], normcilow[inds], normcihi[inds],avg_dates[inds])
names(outdf2) <- c("populations", "scores", "lower", "upper","times")
write.table(outdf2, file=paste(outpath, outtype, ".betase.norm.pval.", pval, sep=""), quote=FALSE, row.names=FALSE)







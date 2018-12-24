# Rscript plot_prs.r 


library("ggplot2")

args<-commandArgs(TRUE)
mpath <- args[1]
pval <- args[2]

fname <- paste("height_pruned.betase.norm.summary.pval.", pval, sep="")

f <- read.table(paste(mpath,fname, sep=""))

giant <- data.frame()
nealelab <- data.frame()

for (i in c(1:nrow(f))){
  giant <- rbind(giant, f[i,1:4])
  nealelab <- rbind(nealelab, f[i,c("V1", "V5", "V6", "V7")])
}

## make plottable data frame segregating by GWAS study
orderx <- c(4,3,1,2) ## reorder populations for creating dataframes

giant <- cbind(giant, "GIANT")
nealelab <- cbind(nealelab, "UK Biobank")

names(giant) <- c("population", "score", "cilo", "cihi", "GWAS")
names(nealelab) <- c("population", "score", "cilo", "cihi", "GWAS")

giant <- giant[orderx,]
nealelab <- nealelab[orderx,]


df <- rbind(giant, nealelab)
names(df) <- c("population", "score", "cilo", "cihi", "GWAS")

## reorder factors for plotting
orderp <- c(1,2,4,3)
df$population <- factor(df$population, levels = df$population[rep(orderp,1)])

### plotting

p <- ggplot(df,aes(x=population,y=score, colour=GWAS)) + 
  geom_pointrange(aes(ymin=cilo, ymax=cihi),position = position_dodge(width = 0.6), size=0.5) + 
  theme_bw() +
  #ylim(-1,1) +
  xlab("") + ylab("Polygenic score") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  theme(panel.border = element_rect(linetype = 5, fill = NA),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text = element_text(size=11, face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=10) )

ggsave(p, file=paste(mpath, "Polygenic_score_2gwas_pval_", pval, ".pdf",sep=""),width=11.5, height=6.5,dpi=300)

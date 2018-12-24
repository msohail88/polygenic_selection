library("ggplot2")

mpath <- ""
fname <- paste("modern")

setwd(mpath)
f <- read.table(paste(fname, sep=""))

prefixes <- c("giant", "nealelab", "poru", "all.nopc","all.pc", "brit.nopc", "brit.pc", "NG2015.sibs", "ukbsibs", "ukbwbsibs")

giant <- data.frame()
nealelab <- data.frame()
poru <- data.frame()
allnopc <- data.frame()
allpc <- data.frame()
britnopc <- data.frame()
britpc <- data.frame()
ng2015.sibs <- data.frame()
ukbsibs <- data.frame()
ukbwbsibs <- data.frame()

for (i in c(1:nrow(f))){
  giant <- rbind(giant, f[i,1:4])
  nealelab <- rbind(nealelab, f[i,c("V1", "V5", "V6", "V7")])
  poru <- rbind(poru, f[i,c("V1", "V8","V9","V10")])
  allnopc <- rbind(allnopc, f[i,c("V1", "V11","V12","V13")])
  allpc <- rbind(allpc, f[i,c("V1", "V14","V15","V16")])
  britnopc <- rbind(britnopc, f[i,c("V1", "V17","V18","V19")])
  britpc <- rbind(britpc, f[i,c("V1", "V20","V21","V22")])
  ng2015.sibs <- rbind(ng2015.sibs, f[i,c("V1", "V23","V24","V25")])
  ukbsibs <- rbind(ukbsibs, f[i,c("V1", "V26","V27","V28")])
  ukbwbsibs <- rbind(ukbwbsibs, f[i,c("V1", "V29","V30","V31")])
}

## make plottable data frame segregating by GWAS study
orderx <- c(4,3,1,2) ## reorder populations for creating dataframes

britpc <- cbind(britpc, "UKB WB 10 PCs")
britnopc <- cbind(britnopc, "UKB WB no PCs")
allpc <- cbind(allpc, "UKB all 10 PCs")
allnopc <- cbind(allnopc, "UKB all no PCs")
giant <- cbind(giant, "GIANT")
poru <- cbind(poru, "UKB Loh")
nealelab <- cbind(nealelab, "UKB Neale")
ng2015.sibs <- cbind(ng2015.sibs, "NG2015 sibs")
ukbsibs <- cbind(ukbsibs, "UKB sibs all")
ukbwbsibs <- cbind(ukbwbsibs,"UKB sibs WB")

names(britpc) <- c("population", "score", "cilo", "cihi", "GWAS")
names(britnopc) <- c("population", "score", "cilo", "cihi", "GWAS")
names(allpc) <- c("population", "score", "cilo", "cihi", "GWAS")
names(allnopc) <- c("population", "score", "cilo", "cihi", "GWAS")
names(giant) <- c("population", "score", "cilo", "cihi", "GWAS")
names(poru) <- c("population", "score", "cilo", "cihi", "GWAS")
names(nealelab) <- c("population", "score", "cilo", "cihi", "GWAS")
names(ng2015.sibs) <- c("population", "score", "cilo", "cihi", "GWAS")
names(ukbsibs) <- c("population", "score", "cilo", "cihi", "GWAS")
names(ukbwbsibs) <- c("population", "score", "cilo", "cihi", "GWAS")

britpc <- britpc[orderx,]
britnopc <- britnopc[orderx,]
allpc <- allpc[orderx,]
allnopc <- allnopc[orderx,]
giant <- giant[orderx,]
poru <- poru[orderx,]
nealelab <- nealelab[orderx,]
ng2015.sibs <- ng2015.sibs[orderx,]
ukbsibs <- ukbsibs[orderx,]
ukbwbsibs <- ukbwbsibs[orderx,]

df <- rbind(giant, nealelab, poru, allnopc, allpc, britnopc, britpc, ng2015.sibs, ukbsibs, ukbwbsibs)
names(df) <- c("population", "score", "cilo", "cihi", "GWAS")

## reorder factors for plotting
orderp <- c(1,2,4,3)
df$population <- factor(df$population, levels = df$population[rep(orderp,10)])

### plotting
p <- ggplot(df,aes(x=population,y=score)) +
  facet_wrap( ~ GWAS, ncol=5) +
  geom_pointrange(aes(ymin=cilo, ymax=cihi),position = position_dodge(width = 0.6), size=0.4) + 
  theme_bw() +
  ylim(-1.5,1.5) +
  xlab("") + ylab("Polygenic height score") +
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

ggsave(p, file=paste("fig1_figsupplement4.pdf",sep=""),width=11.5, height=6.5,dpi=300)




library("ggplot2")
mpath <- ""
fname <- "global.modern"

setwd(mpath)
f <- read.table(paste(fname, sep=""))

giant <- data.frame()
poru <- data.frame()
nealelab <- data.frame()

for (i in c(1:nrow(f))){
  giant <- rbind(giant, f[i,1:4])
  poru <- rbind(poru, f[i,c("V1", "V5", "V6", "V7")])
  nealelab <- rbind(nealelab, f[i,c("V1", "V8","V9","V10")])
}

## make plottable data frame segregating by GWAS study
orderx <- c(12,13,11, 1,2,3,4,9,10,5,6,7,8) ## reorder populations for plotting

giant <- cbind(giant, "GIANT")
poru <- cbind(poru, "UKB Loh")
nealelab <- cbind(nealelab, "UKB Nealelab")

names(giant) <- c("population", "score", "cilo", "cihi", "GWAS")
names(poru) <- c("population", "score", "cilo", "cihi", "GWAS")
names(nealelab) <- c("population", "score", "cilo", "cihi", "GWAS")

giant <- giant[orderx,]
poru <- poru[orderx,]
nealelab <- nealelab[orderx,]

giant$population <- factor(giant$population, levels = giant$population[orderx])
poru$population <- factor(poru$population, levels = poru$population[orderx])
nealelab$population <- factor(nealelab$population, levels = nealelab$population[orderx])

groups <- c("Ancient Europe", "Modern Europe", "East Asia", "South Asia", "Africa")
gr_order <- c(rep(groups[1],3), rep(groups[2],4),rep(groups[3],2),rep(groups[4],2),rep(groups[5],2))

df <- rbind(giant, poru, nealelab)
df <- cbind(df, rep(gr_order,3))
names(df) <- c("population", "score", "cilo", "cihi", "GWAS", "group")

orderg <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
df$group <- factor(df$group, levels = df$group[rep(orderg,3)])

### plotting
p <- ggplot(df,aes(x=population,y=score, colour=GWAS)) + 
  geom_pointrange(aes(ymin=cilo, ymax=cihi), position=position_dodge(width = 0.6), size=0.5, fatten=1.4) + 
  theme_bw() +
  ylim(-4,4) +
  xlab("") + ylab("Polygenic score") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_grid(. ~ group,scales="free",space="free") +
  theme(panel.border = element_rect(linetype = 5, fill = NA),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text = element_text(size=8, face="bold"),
        axis.text=element_text(size=8),
        axis.title=element_text(size=11),
        legend.text=element_text(size=9) )

ggsave(p, file="fig1_figsupplement6.pdf",width=5, height=3.5, dpi=300)
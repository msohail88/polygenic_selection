
library("ggplot2")

mpath <- ""
setwd(mpath)

###### Add modern data ##########
fname <- paste("modern.0.01", sep="")

f <- read.table(paste(fname, sep=""))

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
df$population <- factor(df$population, levels = df$population[rep(orderp,2)])

##
groups <- c("Modern Europe")
gr_order <- c(rep(groups[1],4))

df <- rbind(giant, nealelab)
df <- cbind(df, rep(gr_order,2))
names(df) <- c("population", "score", "cilo", "cihi", "GWAS", "group")

###### Add ancient data ##########
fname <- paste("ancient.0.01", sep="")

f <- read.table(paste(fname, sep=""))

giant <- data.frame()
nealelab <- data.frame()

for (i in c(1:nrow(f))){
  giant <- rbind(giant, f[i,1:4])
  nealelab <- rbind(nealelab, f[i,c("V1", "V5", "V6", "V7")])
}

## make plottable data frame segregating by GWAS study
orderx <- c(1,2,3) ## reorder populations for creating dataframes

giant <- cbind(giant, "GIANT")
nealelab <- cbind(nealelab, "UK Biobank")

names(giant) <- c("population", "score", "cilo", "cihi", "GWAS")
names(nealelab) <- c("population", "score", "cilo", "cihi", "GWAS")

giant <- giant[orderx,]
nealelab <- nealelab[orderx,]


df2 <- rbind(giant, nealelab)
names(df2) <- c("population", "score", "cilo", "cihi", "GWAS")

## reorder factors for plotting
orderp <- c(2,3,1)
df2$population <- factor(df2$population, levels = df2$population[rep(orderp,2)])

##
groups <- c("Ancient Europe")
gr_order <- c(rep(groups[1],3))

df2 <- rbind(giant, nealelab)
df2 <- cbind(df2, rep(gr_order,2))
names(df2) <- c("population", "score", "cilo", "cihi", "GWAS", "group")

###### combine modern and ancient DNA #######
df_final <- rbind(df, df2)

## reorder factors for plotting
orderp <- c(1,2,4,3,5,6,8,7,11,9,10,14,12,13)
df_final$population <- factor(df_final$population, levels = df_final$population[orderp])

### plotting

p <- ggplot(df_final,aes(x=population,y=score, colour=GWAS)) + 
  geom_pointrange(aes(ymin=cilo, ymax=cihi),position = position_dodge(width = 0.6), size=0.8) + 
  theme_bw() +
  ylim(-4,4) +
  xlab("") + ylab("Polygenic height score") +
  geom_hline(yintercept=0,linetype="dashed") +
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid")) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_grid(. ~ group,scales="free",space="free") +
  theme(panel.border = element_rect(linetype = 5, fill = NA),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0.15, "lines"),
        strip.text = element_text(size=11, face="bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=13),
        legend.text=element_text(size=10) )

ggsave(p, file=paste("fig1_figsupplement2a.pdf",sep=""),width=9, height=6, dpi=300)
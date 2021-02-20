library(ggplot2)
library(reshape2)
library(scales)

OsHigh <- read.csv("Osativa_high_gc_metaplot.tsv",header=T,sep="\t")
OsLow <- read.csv("Osativa_low_gc_metaplot.tsv",header=T,sep="\t")
SbHigh <- read.csv("Sbicolor_high_gc_metaplot.tsv",header=T,sep="\t")
SbLow <- read.csv("Sbicolor_Low_gc_metaplot.tsv",header=T,sep="\t")

OsHigh <- melt(OsHigh[,2:4])
OsLow <- melt(OsLow[,2:4])
SbHigh <- melt(SbHigh[,2:4])
SbLow <- melt(SbLow[,2:4])

OsHigh$Window <- c(1:60)
OsLow$Window <- c(1:60)
SbHigh$Window <- c(1:60)
SbLow$Window <- c(1:60)

OsHigh$Syntelogs <- c('High GC')
OsLow$Syntelogs <- c('Low GC')
SbHigh$Syntelogs <- c('High GC')
SbLow$Syntelogs <- c('Low GC')

Os <- rbind(OsHigh,OsLow)
Sb <- rbind(SbHigh,SbLow)

Os$variable <- gsub("_Weighted_mC","",Os$variable)
Sb$variable <- gsub("_Weighted_mC","",Sb$variable)

p <- ggplot(Os) + geom_line(aes(x=Window,y=value,color=variable,
                           linetype=Syntelogs),size=0.8) + 
  theme_bw() + scale_y_continuous("Weighted Methylation",labels=percent,
                                  expand=c(0,0)) +
  scale_x_continuous("Coding Sequences",expand=c(0,0),breaks=c(1,20,40,60),
                     labels=c("-2000bp","Start","Stop","+2000bp")) +
  geom_vline(aes(xintercept=20),linetype="longdash",color="grey65") + 
  geom_vline(aes(xintercept=40),linetype="longdash",color="grey65") +
  scale_color_manual("mC Type",
                     values=c("violetred2","dodgerblue3","cyan3")) +
  ggtitle("O. sativa")
ggsave("Osativa_GC_syntelog_metaplot.pdf",p,device="pdf")

p <- ggplot(Sb) + geom_line(aes(x=Window,y=value,color=variable,
                                linetype=Syntelogs),size=0.8) + 
  theme_bw() + scale_y_continuous("Weighted Methylation",labels=percent,
                                  expand=c(0,0)) +
  scale_x_continuous("Coding Sequences",expand=c(0,0),breaks=c(1,20,40,60),
                     labels=c("-2000bp","Start","Stop","+2000bp")) +
  geom_vline(aes(xintercept=20),linetype="longdash",color="grey65") + 
  geom_vline(aes(xintercept=40),linetype="longdash",color="grey65") +
  scale_color_manual("mC Type",
                     values=c("violetred2","dodgerblue3","cyan3")) +
  ggtitle("S. bicolor")
ggsave("Sbicolor_GC_syntelog_metaplot.pdf",p,device="pdf")

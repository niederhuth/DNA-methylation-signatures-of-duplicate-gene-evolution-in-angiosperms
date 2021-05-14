##############################
#   plots from methylome data (928 Arabidopsis accesions) 
##############################

library(ggplot2)
library(Cairo)

df1 <- read.table("Athaliana/methylpy/results/Athaliana_classified_genes.tsv",header=T,
	sep="\t")[c("Feature","Classification")]
df2 <- read.table("Athaliana/dupgen/results-unique/classified_genes.tsv",header=T,sep="\t")
df3 <- merge(df1,df2,by="Feature")

#Read in the table of 928 accessions
df4 <- read.table("Athaliana/methylpy/results/At_variation.tsv",header=TRUE,sep="\t")
#Convert "NAs" to "Missing"
df4[is.na(df4)] <- "Missing"
#Create an empty dataframe to summarize methylation frequencies
df5 <- data.frame(Feature=character(),gbM=numeric(),teM=numeric(),unM=numeric(),
	Unclassified=numeric(),Missing=numeric())
#Loop over every gene in df1 to build frequency table
for(i in c(1:nrow(df1))){df1[i,]
	#create a temporary dataframe after tabling the classifications in df4
	x <- as.data.frame(table(t(df4[i,])))
	#add these values to each column in df5, if that value is not present, set to "NA"
	df5 <- rbind(df5,data.frame(Feature=df4[i,]$Feature,
		gbM=ifelse(x[x[1]=="gbM",]$Freq,x[x[1]=="gbM",]$Freq,NA)[1],
		teM=ifelse(x[x[1]=="teM",]$Freq,x[x[1]=="teM",]$Freq,NA)[1],
		unM=ifelse(x[x[1]=="unM",]$Freq,x[x[1]=="unM",]$Freq,NA)[1],
		Unclassified=ifelse(x[x[1]=="Unclassified",]$Freq,x[x[1]=="Unclassified",]$Freq,NA)[1],
		Missing=ifelse(x[x[1]=="Missing",]$Freq,x[x[1]=="Missing",]$Freq,NA)[1]))
}
#Sum up the rows. We do not count data that is missing
df5$Total <- rowSums(df5[2:5],na.rm=TRUE)
#Get percentage of each classification
df5$p_gbM <- df5$gbM/df5$Total
df5$p_teM <- df5$teM/df5$Total
df5$p_unM <- df5$unM/df5$Total
df5$p_Unclassified <- df5$Unclassified/df5$Total
#Merge dataframes
df6 <- merge(df3,df5,by="Feature")
#Group genes based on frequency of gbM
df6$gbMgroup <- NA
df6$gbMgroup <- ifelse(df6$p_gbM > 0.75, 5, df6$gbMgroup)
df6$gbMgroup <- ifelse(df6$p_gbM < 0.75, 4, df6$gbMgroup)
df6$gbMgroup <- ifelse(df6$p_gbM < 0.5, 3, df6$gbMgroup)
df6$gbMgroup <- ifelse(df6$p_gbM < 0.25, 2, df6$gbMgroup)
df6$gbMgroup <- ifelse(is.na(df6$p_gbM), 1, df6$gbMgroup)
#Group genes based on frequency of teM
df6$teMgroup <- NA
df6$teMgroup <- ifelse(df6$p_teM > 0.75, 5, df6$teMgroup)
df6$teMgroup <- ifelse(df6$p_teM < 0.75, 4, df6$teMgroup)
df6$teMgroup <- ifelse(df6$p_teM < 0.5, 3, df6$teMgroup)
df6$teMgroup <- ifelse(df6$p_teM < 0.25, 2, df6$teMgroup)
df6$teMgroup <- ifelse(is.na(df6$p_teM), 1, df6$teMgroup)
#Group genes based on frequency of unM
df6$unMgroup <- NA
df6$unMgroup <- ifelse(df6$p_unM > 0.75, 5, df6$unMgroup)
df6$unMgroup <- ifelse(df6$p_unM < 0.75, 4, df6$unMgroup)
df6$unMgroup <- ifelse(df6$p_unM < 0.5, 3, df6$unMgroup)
df6$unMgroup <- ifelse(df6$p_unM < 0.25, 2, df6$unMgroup)
df6$unMgroup <- ifelse(is.na(df6$p_unM), 1, df6$unMgroup)
#Read in KaKs values from KaKs.R script
df7 <- read.csv("../figures_tables/Athaliana/Athaliana_kaks_values.csv",header=TRUE)
#Merge the dataframes
#Some ugly code, but it works
df8 <- data.frame()
for(a in 1:nrow(df6)){
	tmp <- df7[df7$Duplication==df6[a,]$Duplication & df7$Duplicate.1==df6[a,]$Feature | 
			df7$Duplication==df6[a,]$Duplication & df7$Duplicate.2==df6[a,]$Feature, ]
	if(nrow(tmp) == 1){
		df8 <- rbind(df8,cbind(df6[a,],tmp))
	} else if(nrow(tmp) > 1){
		tmp2 <- read.table(paste("Athaliana/dupgen/results-unique/Athaliana.",df6[a,]$Duplication,
								".pairs-unique",sep=""),header=TRUE,sep="\t")
		tmp3 <- read.table(paste("Athaliana/dupgen/results/Athaliana.",df6[a,]$Duplication,
								".pairs",sep=""),header=TRUE,sep="\t")
		tmp4 <- unique(rbind(merge(tmp,tmp2,by.x=c("Duplicate.1","Duplicate.2"),by.y=c("Duplicate.1","Duplicate.2")),
				merge(tmp,tmp2,by.x=c("Duplicate.1","Duplicate.2"),by.y=c("Duplicate.2","Duplicate.1")),
				merge(tmp,tmp3,by.x=c("Duplicate.1","Duplicate.2"),by.y=c("Duplicate.1","Duplicate.2")),
				merge(tmp,tmp3,by.x=c("Duplicate.1","Duplicate.2"),by.y=c("Duplicate.2","Duplicate.1"))))
		df8 <- rbind(df8,cbind(df6[a,],tmp4[tmp4$E.value==min(tmp4$E.value),c(1:11)]))
	}
}
#Clean up the table 
df9 <- df8
colnames(df9) <- c("Feature","Classification","Duplication","gbM","teM","unM","Unclassified","Missing","Total",
					"p_gbM","p_teM","p_unM","p_Unclassified","gbMgroup" ,"teMgroup","unMgroup","Duplicate.1",
					"Duplicate.2","Duplicate.1_Methylation","Duplicate.2_Methylation","Duplication2" ,"Similarity",
					"Dup_Met_Class","Ka","Ks","Ka.Ks","P.Value")

#Labels for plots
labels <- c("0%","<25%","25-50%","50-75%",">75%")
#Output directory
path1 <- "../figures_tables/Athaliana/Diversity/"
if(!file.exists(path1)){
	dir.create(path1)
}
#Make plots
setcolors <- c("#F8766D", "#C49A00","#53B400", "#00B6EB", "#FB61D7")

p <- ggplot(df9[!is.na(df9$gbMgroup),]) + 
	geom_density(aes(x=Ks,color=as.character(gbMgroup)),size=0.5) + 
	scale_color_discrete(labels=labels) +
	theme_bw()+
	scale_x_continuous(expand=c(0.05,0.05), limits = c(0,5)) + 
	scale_y_continuous(expand=c(0,0), limits=c(0,0.7))
ggsave(paste(path1,"gbM_Ks_density.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$teMgroup),]) + 
	geom_density(aes(x=Ks,color=as.character(teMgroup)),size=0.5) + 
	scale_color_discrete(labels=labels) +
	theme_bw()+
	scale_x_continuous(expand=c(0.05,0.05), limits = c(0,5)) + 
	cale_y_continuous(expand=c(0,0), limits=c(0,0.7))
ggsave(paste(path1,"teM_Ks_density.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$unMgroup),]) + 
	geom_density(aes(x=Ks,color=as.character(unMgroup)),size=0.5) + 
	scale_color_discrete(labels=labels)+
	theme_bw()+
	scale_x_continuous(expand=c(0.05,0.05), limits = c(0,5)) + 
	scale_y_continuous(expand=c(0,0), limits=c(0,0.7))
ggsave(paste(path1,"unM_Ks_density.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$gbMgroup),]) + 
	geom_boxplot(aes(y=Ks,x=as.character(gbMgroup)),col=setcolors,size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + 
	scale_y_continuous(expand=c(0,0), limits = c(0,5))
ggsave(paste(path1,"gbM_Ks_boxplot.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$teMgroup),]) + 
	geom_boxplot(aes(y=Ks,x=as.character(teMgroup)),col=setcolors,size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + 
	scale_y_continuous(expand=c(0,0), limits = c(0,5))
ggsave(paste(path1,"teM_Ks_boxplot.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$unMgroup),]) + 
	geom_boxplot(aes(y=Ks,x=as.character(unMgroup)),col=setcolors,size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + 
	scale_y_continuous(expand=c(0,0), limits = c(0,5))
ggsave(paste(path1,"unM_Ks_boxplot.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$gbMgroup),]) + 
	geom_boxplot(aes(y=Ka.Ks,x=as.character(gbMgroup)),col=setcolors,size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + 
	scale_y_continuous(expand=c(0,0),limits=c(0,2))
ggsave(paste(path1,"gbM_KaKs_boxplot.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$teMgroup),]) + 
	geom_boxplot(aes(y=Ka.Ks,x=as.character(teMgroup)),col=setcolors,size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + 
	scale_y_continuous(expand=c(0,0),limits=c(0,2))
ggsave(paste(path1,"teM_KaKs_boxplot.pdf",sep=""),p)

p <- ggplot(df9[!is.na(df9$unMgroup),]) + 
	geom_boxplot(aes(y=Ka.Ks,x=as.character(unMgroup)),col=setcolors,size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + 
	scale_y_continuous(expand=c(0,0),limits=c(0,2))
ggsave(paste(path1,"unM_KaKs_boxplot.pdf",sep=""),p)

df9 <- subset(df9, !Duplication=="singletons")
p <- ggplot(df9) + geom_violin(aes(y=Ks, x=as.character(teMgroup)), size=1) +
  scale_x_discrete(labels = labels) + 
  theme_bw() + 
  scale_y_continuous(expand=c(0.1,0.1), limits = c(0,5))+
  facet_wrap(~df9$Duplication, nrow = 5)
ggsave(paste(path1,"Ks_violin.pdf",sep=""),p)


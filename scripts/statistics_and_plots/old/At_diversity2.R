
##############################
#   plots from methylome data (928 Arabidopsis accesions) 
##############################

library(ggplot2)
library(Cairo)

#setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/data/Athaliana")

df <- read.table("methylpy/results/Athaliana_classified_genes.tsv",header=T,
	sep="\t")[c("Feature","Classification")]
df2 <- read.table("dupgen/results-unique/classified_genes.tsv",header=T,sep="\t")
df3 <- merge(df,df2,by="Feature")

#Read in the table of 928 accessions
df4 <- read.table("methylpy/results/At_variation.tsv",header=TRUE,sep="\t")
#Convert "NAs" to "Missing"
df4[is.na(df4)] <- "Missing"
#Create an empty dataframe to summarize methylation frequencies
df5 <- data.frame(Feature=character(),gbM=numeric(),teM=numeric(),unM=numeric(),
	Unclassified=numeric(),Missing=numeric())
#Loop over every gene in df
for(i in c(1:nrow(df))){df[i,]
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
df7 <- read.csv("../../figures_tables/Athaliana/Athaliana_KaKs_values.csv",header=TRUE,
	row.names=1)[c(1,5,6,7,8)]
#Merge files
df8 <- merge(df6,df7,by="Feature",all=TRUE)
#Create new df of only duplicate pairs
df9 <- df8[df8$Duplication != "singelton" & df8$Duplication != "unclassified",]
#Labels for plots
labels <- c("0%","<25%","25-50%","50-75%",">75%")

ggplot(df9[!is.na(df9$gbMgroup),]) + geom_density(aes(x=Ks,color=as.character(gbMgroup)),size=0.5) + 
	scale_color_discrete(labels=labels) +
	theme_bw()+
	scale_x_continuous(expand=c(0.05,0.05)) + scale_y_continuous(expand=c(0,0))

ggplot(df9[!is.na(df9$teMgroup),]) + geom_density(aes(x=Ks,color=as.character(teMgroup)),size=0.5) + 
	scale_color_discrete(labels=labels) +
	theme_bw()+
	scale_x_continuous(expand=c(0.05,0.05)) + scale_y_continuous(expand=c(0,0))

ggplot(df9[!is.na(df9$unMgroup),]) + geom_density(aes(x=Ks,color=as.character(unMgroup)),size=0.5) + 
	scale_color_discrete(labels=labels)+
	theme_bw()+
	scale_x_continuous(expand=c(0.05,0.05)) + scale_y_continuous(expand=c(0,0))

ggplot(df9[!is.na(df9$gbMgroup),]) + geom_boxplot(aes(y=Ks,x=as.character(gbMgroup)),size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + scale_y_continuous(expand=c(0,0))

ggplot(df9[!is.na(df9$teMgroup),]) + geom_boxplot(aes(y=Ks,x=as.character(teMgroup)),size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + scale_y_continuous(expand=c(0,0))

ggplot(df9[!is.na(df9$unMgroup),]) + geom_boxplot(aes(y=Ks,x=as.character(unMgroup)),size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + scale_y_continuous(expand=c(0,0))

ggplot(df9[!is.na(df9$gbMgroup),]) + geom_boxplot(aes(y=Ka.Ks,x=as.character(gbMgroup)),size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + scale_y_continuous(expand=c(0,0),limits=c(0,2))

ggplot(df9[!is.na(df9$teMgroup),]) + geom_boxplot(aes(y=Ka.Ks,x=as.character(teMgroup)),size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + scale_y_continuous(expand=c(0,0),limits=c(0,2))

ggplot(df9[!is.na(df9$unMgroup),]) + geom_boxplot(aes(y=Ka.Ks,x=as.character(unMgroup)),size=1,notch=TRUE) + 
	theme_bw()+
	theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
	theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
	scale_x_discrete(labels=labels) + scale_y_continuous(expand=c(0,0),limits=c(0,2))
library(DESeq2)


#Function to calculate tau for expression specificity. See Yanai et al. 2015 Bioinformatics
tauSpecificity <- function(x){
	i = 0
	for(j in c(1:length(x))){
		i = i + 1 - (log2(x[j])/log2(max(x)))
	}
	y = i/(length(x) - 1)
	return(y)
}

setwd("C:/Users/kench/Dropbox/ANALYSIS/GeneDuplication_V2/data") #new laptop
for(i in c("Athaliana")){
	#Input data path
	path1=paste(i,"/rna",sep="")
	#Output data path
	path2=paste("../figures_tables/",i,sep="")
	#Get list of samples from directories in rna
	samples <- list.files(path=path1)
	files <- list.files(path=path1,pattern="_ReadsPerGene.out.tab",recursive=TRUE)
	#Get gene names
	# genes <- read.table(paste(path1,files[1],sep="/"),header=FALSE,sep="\t",skip=4)[,1]
	# #Modify gene names
	# genes <- gsub(".Araport11.447","",genes)
	# #Read in files and build the table
	# df <- do.call(cbind,lapply(files,
	# 	function(fn)read.table(paste(path1,fn,sep="/"),header=FALSE, sep="\t", skip=4)[,2]))
	# #Add genes as rownames
	# row.names(df)<- genes
	# #Add samples as column names 
	# colnames(df) <- samples
	# #Output as a table
	# #write.csv(df,paste(path2,"/",i,"_gene_exp.csv",sep=""),quote=FALSE,row.names=TRUE)
	# #Create dataframe describing samples & replicates
	# coldata <- data.frame(row.names=colnames(df),condition=gsub("-.$","",samples))
	# #Create a DESeqDataSet from the matrix and coldata
	# dds <- DESeqDataSetFromMatrix(countData=df, colData=coldata, design = ~ condition)
	# #Normalize data
	# dds <- estimateSizeFactors(dds)
	# #Extract out normalized counts
	# df2 <- data.frame(counts(dds,normalized=TRUE), check.names = F)
	# #write.csv(df2,paste(path2,"/",i,"_gene_exp_normalized.csv",sep=""),quote=FALSE,row.names=TRUE)
	# #Create an empty dataframe to get means of different conditions
	# df3 <- data.frame(row.names=row.names(df2))
	# #Create a new table, where we get the mean of the conditions with replicates
	# #Since some conditions only had one replicate, we just keep the original value for those conditions
	# for(x in unique(coldata$condition)){
	# 	if(ncol(data.frame(df2[,grepl(x, colnames(df2))])) > 1){
	# 		df3[,x] <- rowMeans(df2[,grepl(x, colnames(df2))])
	# 	} else {
	# 		df3[,x] <- df2[,grepl(x, colnames(df2))]
	# 	}
	# }
	#We need to add some number to the data, so that we can log transform samples with 0 read counts
	#without getting infinite values
	
	df3 <- read.csv("./Athaliana/Athaliana_normalized_expression.csv", header =TRUE, check.names = F)
	rownames(df3) <- (df3[,1])
	df3 <- df3[,-1]
	df3 <- df3 + 1
	#We now create a new data frame, applying the tau tissue specificity 
	df4 <- data.frame(tau=apply(df3,1,tauSpecificity))
	#Add gene names as a column called features
	df4$Feature <- row.names(df4)
	write.csv(df4,paste(path2,"/",i,"_gene_Tau_05-10.csv",sep=""),quote=FALSE,row.names=TRUE)
	#Read in KaKs table, with duplicate gene pairs and other relevant info
	df5 <- read.csv(paste("../figures_tables/",i,"/",i,"_KaKs_values.csv",sep=""),
		header=TRUE,row.names=1)
	#Merge the data frames together.
	df6 <- merge(df5,df4,by="Feature")
	write.csv(df6,paste(path2,"/",i,"_gene_KaKs_Tau_05-10.csv",sep=""),quote=FALSE,row.names=TRUE)
}

#### Making figures for tau ####

library(ggplot2)
library(plyr)

setwd("C:/Users/kench/Dropbox/ANALYSIS/GeneDuplication_V2/data")
df7 <- read.table("./Athaliana/methylpy/results/Athaliana_classified_genes.tsv",sep="\t", header=TRUE) #Methylation classification
df8 <- merge(df4,df7[,c(1,30)],by.x="Feature",by.y="Feature") #combining methylation classification and Tau
df8 <- df8[df8$Classification != 'Unclassified', ]
df8 <- df8[complete.cases(df8),] # removes NA data
write.csv(df8, "../figures_tables/Athaliana/Athaliana_genes_Tau_05-10.csv")
#df9 <- read.csv("../figures_tables/Athaliana/Athaliana_Duplicate_pair_met2.csv", header=TRUE, check.names=F) #Duplicate similarity
df9 <- read.csv("../figures_tables/Athaliana/DuplicateSimilarity_sorted.csv", header=TRUE, check.names=F) # manually made file without undetermined 
#df10 <- merge(df8, df9[,c(2,7)],by.x="Feature",by.y="Duplicate.1")
df11 <- merge(df9, df8[,c(1:2)], by.x="Duplicate.1",by.y="Feature")
ordered <- df11[order(df11$id), ]

df12 <- merge(df9, df8[,c(1:2)], by.x="Duplicate.2",by.y="Feature")
ordered2 <- df12[order(df12$id), ]

df13 <- merge(ordered[,c(1,2,3,6,8,9)],ordered2[,c(1,2,5, 8,9)], by = "id") #has Tau for dup1 and Dup2 and added duplication type from ordered
write.csv(df13, "../figures_tables/Athaliana/DuplicateSimilarity_Tau_SeperateDuplicates.csv")
df14 <- df13[,c(1,4,5,6, 10)]
df14$absdiff <- abs(df14$tau.x-df14$tau.y)

write.csv(df14, "DuplicateSimilarity_Tau_05-10.csv")

#### Athaliana RNA expression figures ####
#This creates a facet graph for each classification combination but gives a wrong sense of directionality
library(GGally) # to use ggparcoord for parallel coordinated chart
#columns should be the values on which parallel plottign is to be performed
#Group column is the x-axis
#scale should be set to globalminmax for no transformation
p <- ggparcoord(df14, scale ="globalminmax", boxplot = TRUE, columns=4:5, groupColumn = 3 ) + 
  facet_wrap(~Classification.x)
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_ggally.pdf" ,p,device="pdf", width = 6, height = 3)

#This creates the gghalves plot with boxplot + jitter showing a good representation of absolute difference in tau
library(gghalves)

df14$Duplication <- factor(df14$Duplication, levels = c("wgd", "tandem", "proximal", "transposed", "dispersed"))
p <- ggplot() +
  theme_bw() +
  geom_half_boxplot(
    data = df14, 
    aes(x = Classification.x, y = absdiff, fill = Classification.x), side = "l", errorbar.draw = TRUE,
    outlier.color = NA) +
  geom_half_point(
    data = df14, 
    aes(x = Classification.x, y = absdiff, fill = Classification.x, color = Classification.x), side = "r", size =0.5) 

ggsave("../figures_tables/Athaliana/Ath_Tau_AbsDiff_gghalves.pdf" ,p,device="pdf", width = 6, height = 3)

p <- ggplot() +
  theme_bw() +
  geom_half_boxplot(
    data = df8, 
    aes(x = Classification, y = tau, fill = Classification), side = "l", errorbar.draw = TRUE,
    outlier.color = NA) +
  geom_half_point(
    data = df8, 
    aes(x = Classification, y = tau, fill = Classification, color = Classification), side = "r", size =0.5) 
ggsave("../figures_tables/Athaliana/Ath_Tau_metclass_gghalves.pdf" ,p,device="pdf", width = 4, height = 3)

#Plot for gbM genes that do or do not switch
#For this I make three different input files that has all gbM genes including undetermined, gbm-gbm, gbm-tem, and gbm-unm.

df20 <- read.csv("../figures_tables/Athaliana/AthRNA_switching_gbm.csv", header =TRUE) # manually made file used as input
df20$Classification.x <- factor(df20$Classification.x, levels = c("All gbM", "gbM-gbM", "gbM-teM", "gbM-unM")) 
p <- ggplot(df20, aes(x=Classification.x, y=tau.x, fill = Classification.x)) + geom_boxplot() +
  #geom_point(alpha=0.5, colour="black", shape=21) +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, size=10), axis.text.y=element_text(color="black", size=10), 
        axis.title.y=element_text(size=12), legend.text=element_text(size=10), panel.border=element_rect(colour="#999999", fill=NA), 
        panel.background=element_rect(fill="white"), panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
        strip.background=element_rect(fill="white"), strip.text.y=element_text(angle=0, hjust=0, size=10),
        plot.title=element_text(size=12)) + 
  scale_fill_brewer(palette="Paired") + 
  labs(x="", y="Tau - gbM genes  \n")
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_gbM.pdf" ,p,device="pdf", width = 4, height = 3)

q <- ggplot() +
  theme_bw() +
  geom_half_boxplot(
    data = df20, 
    aes(x = Classification.x, y = tau.x, fill = Classification.x), side = "l", errorbar.draw = TRUE,
    outlier.color = NA) +
  geom_half_point(
    data = df20, 
    aes(x = Classification.x, y = tau.x, fill = Classification.x, color = Classification.x), side = "r", size =0.5) 
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_gbM_gghalves.pdf" ,q,device="pdf", width = 4, height = 3)

#Plot for teM genes that do or donot switch
df20 <- read.csv("../figures_tables/Athaliana/AthRNA_switching_tem.csv", header =TRUE) # manually made file used as input
df20$Classification.x <- factor(df20$Classification.x, levels = c("All teM", "teM-teM", "teM-gbM", "teM-unM")) 
p <- ggplot(df20, aes(x=Classification.x, y=tau.x, fill = Classification.x)) + geom_boxplot() +
  #geom_point(alpha=0.5, colour="black", shape=21) +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, size=10), axis.text.y=element_text(color="black", size=10), 
        axis.title.y=element_text(size=12), legend.text=element_text(size=10), panel.border=element_rect(colour="#999999", fill=NA), 
        panel.background=element_rect(fill="white"), panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
        strip.background=element_rect(fill="white"), strip.text.y=element_text(angle=0, hjust=0, size=10),
        plot.title=element_text(size=12)) +
  scale_fill_brewer(palette="Paired") + 
  labs(x="", y="Tau - teM genes  \n")
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_teM.pdf" ,p,device="pdf", width = 4, height = 3)

q <- ggplot() +
  theme_bw() +
  geom_half_boxplot(
    data = df20, 
    aes(x = Classification.x, y = tau.x, fill = Classification.x), side = "l", errorbar.draw = TRUE,
    outlier.color = NA) +
  geom_half_point(
    data = df20, 
    aes(x = Classification.x, y = tau.x, fill = Classification.x, color = Classification.x), side = "r", size =0.5) 
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_teM_gghalves.pdf" ,q,device="pdf", width = 4, height = 3)

#Plot for unM genes that do or donot switch
df20 <- read.csv("../figures_tables/Athaliana/AthRNA_switching_unm.csv", header =TRUE) # manually made file used as input
df20$Classification.x <- factor(df20$Classification.x, levels = c("All unM", "unM-unM", "unM-gbM", "unM-teM")) 
p <- ggplot(df20, aes(x=Classification.x, y=tau.x, fill = Classification.x)) + geom_boxplot() +
  #geom_point(alpha=0.5, colour="black", shape=21) +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, size=10), axis.text.y=element_text(color="black", size=10), 
        axis.title.y=element_text(size=12), legend.text=element_text(size=10), panel.border=element_rect(colour="#999999", fill=NA), 
        panel.background=element_rect(fill="white"), panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
        strip.background=element_rect(fill="white"), strip.text.y=element_text(angle=0, hjust=0, size=10),
        plot.title=element_text(size=12)) +
  scale_fill_brewer(palette="Paired") + 
  labs(x="", y="Tau - unM genes  \n")
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_unM.pdf" ,p,device="pdf", width = 4, height = 3)

q <- ggplot() +
  theme_bw() +
  geom_half_boxplot(
    data = df20, 
    aes(x = Classification.x, y = tau.x, fill = Classification.x), side = "l", errorbar.draw = TRUE,
    outlier.color = NA) +
  geom_half_point(
    data = df20, 
    aes(x = Classification.x, y = tau.x, fill = Classification.x, color = Classification.x), side = "r", size =0.5)

ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_unM_gghalves.pdf" ,q,device="pdf", width = 4, height = 3)

df40 <- read.table("Athaliana/dupgen/results-unique/classified_genes.tsv", sep = "\t", header = T)
df41 <- merge (df8, df40, by.x="Feature",by.y="Feature")
df41 <- df41[df41$Duplication != "unclassified",]
df41 <- df41[df41$Duplication != "singletons",]


df41$Duplication <- factor(df41$Duplication, levels = c("wgd", "tandem", "proximal", "transposed", "dispersed")) 
p <- ggplot(df41, aes(x=Duplication, y=tau, fill=Duplication)) + geom_boxplot() +
  #geom_point(alpha=0.5, colour="black", shape=21) +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, size=10), axis.text.y=element_text(color="black", size=10), 
        axis.title.y=element_text(size=12), legend.text=element_text(size=10), panel.border=element_rect(colour="#999999", fill=NA), 
        panel.background=element_rect(fill="white"), panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
        strip.background=element_rect(fill="white"), strip.text.y=element_text(angle=0, hjust=0, size=10),
        plot.title=element_text(size=12)) + 
  labs(x="", y="Tau \n") +
  scale_fill_brewer(palette="Paired") 
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_Dupclass.pdf" ,p,device="pdf", width = 4, height = 3) 

q <- ggplot(df41, aes(x=Duplication, y=tau, fill=Duplication)) + geom_boxplot() +
  #geom_point(alpha=0.5, colour="black", shape=21) +
  theme(axis.text.x=element_text(color="black", angle=45, hjust=1, size=10), axis.text.y=element_text(color="black", size=10), 
        axis.title.y=element_text(size=12), legend.text=element_text(size=10), panel.border=element_rect(colour="#999999", fill=NA), 
        panel.background=element_rect(fill="white"), panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
        strip.background=element_rect(fill="white"), strip.text.y=element_text(angle=0, hjust=0, size=10),
        plot.title=element_text(size=12)) + 
  labs(x="", y="Tau \n") +
  scale_fill_brewer(palette="Paired") +
  facet_grid(~Classification)
ggsave("../figures_tables/Athaliana/Ath_Tau_Switching_Dupclass_metclass.pdf" ,q,device="pdf", width = 6, height = 3)


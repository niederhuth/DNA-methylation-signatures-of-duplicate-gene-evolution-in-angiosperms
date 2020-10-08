library(ggplot2)
library(reshape2)
library(scales)

#List of species
species=c("Athaliana")

#Read in orthogroup gene counts for all species & reformat
path1=paste("orthofinder/orthofinder/",dir("orthofinder/orthofinder/"),sep="")
geneCounts <- read.table(paste(path1,"/Orthogroups/Orthogroups.GeneCount.tsv",sep=""),
	header=T,sep="\t",row.names=1)
#Create two new temporary dataframes to identify missing orthogroups and single-copy orthogroups
df3 <- df2 <- subset(geneCounts,select = -Total )
#If a species is missing an orthogroup label it as a 0, if the orthogroup is present
#label it as a 1, regardless of how many genes are in that orthogroup
for(i in colnames(df2)){
  df2[i] <- ifelse(geneCounts[i]==0,0,1)
}
#Sum up the row to get how many species have or are missing an orthogroup 
df2$Total <- rowSums(df2)
#If an orthogroup is single-copy in a species, label it as a 1, otherwise label it a 0.
for(i in colnames(df3)){
  df3[i] <- ifelse(geneCounts[i]==1,1,0)
}
#Get the percentage of species where an orthogroup is single-copy
df3$Total <- rowSums(df3)/ncol(df3)
#Plot out the distribution of species per orthogroup. We will use this to decide what the
#the minimal number of species at which we classify an orthogroup as being part of the core
#angiosperm gene set. This is similar to what was done in Li et al. 2016
p <- ggplot(df2) + 
	geom_histogram(aes(x=Total,fill=cut(Total,c(50,58))),bins=58) + 
	theme_bw() + theme(legend.position="none") + 
	scale_y_continuous("Number of Orthogroups",expand=c(0,0)) + 
	scale_x_continuous("Number of Species in Orthogroup",expand=c(0,0))
#Make a dataframe of "core" angiosperm genes
coreGenes <- df3[row.names(df3) %in% row.names(df2[df2$Total >= 51,]),]
#Make a dataframe of single-copy core angiosperm genes. We will allow for a certain percentage of these
#to be present as multicopy in some species, due to recent WGD and other duplication events
singleCopy <- geneCounts[row.names(geneCounts) %in%  row.names(coreGenes[coreGenes$Total >= 0.75,]),]

#Iterate over each species
for(a in species){
	#Read in the list of genes classified by methylation status
	df1 <- read.table(paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep=""),sep="\t",header=TRUE)
	#Set classification to character, so that R doesn't screw it up
	df1$Classification <- as.character(df1$Classification)
	#Change "NA" to "Missing"
	df1$Classification <- ifelse(is.na(df1$Classification),"Missing",df1$Classification)
	#Read in the orthogroups for that species
	df2 <- read.table(paste(a,"/ref/mcscanx/",a,"_orthogroups.tsv",sep=""),header=FALSE,sep="\t")
	#Merge DNA methylation classification & orthogroups
	df3 <- merge(df2,df1[c(1,23)],by.x="V1",by.y="Feature")
	#Table the values so that we can get the count of each classification for each orthogroup
	df4 <- data.frame(table(df3[c(2,3)]))
	#Build a new table with each column being the DNA methylation classification
	#and each row an orthogroup
	df5 <- df4[df4$Classification=="gbM",c(1,3)]
	df5 <- merge(df5,df4[df4$Classification=="TE-like",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Unmethylated",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Unclassified",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Missing",c(1,3)],by="V2")
	#Rename columns
	colnames(df5) <- c("Orthogroup","gbM","TE-like","Unmethylated","Unclassified","Missing")
	#Sum up each row
	df5$Total <- rowSums(df5[c(2:6)])
	#Reformat for easier plotting
	df6 <- melt(df5)
}




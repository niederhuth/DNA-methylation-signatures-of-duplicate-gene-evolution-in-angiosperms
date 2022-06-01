#Script for normalizing count data and combining replicates

library(DESeq2)
library(reshape2)

species=c("Athaliana","Gmax","Pvulgaris","Sbicolor")

for(i in species){
	#Input data path
	path1=paste(i,"/rna",sep="")
	#Output data path
	path2=paste("../figures_tables/",i,sep="")
	#Check if Gene_Expression.tsv file exists -> preprocessed JGI datasets
	#Output a reformatted table
	if(!file.exists(path1=paste(i,"/Gene_Expression.tsv"))){
		df1 <- read.table(path1,header=F,sep="\t")
		df1 <- df1[c(1,3,2)]
		colnames(df1) <- c("Feature","variable","value")
		df2 <- dcast(df1,formula=Feature~variable,fun.aggregate=sum,value.var="value")
		row.names(df2) <- df2$Feature
		df3 <- df2[2:ncol(df2)]
	#Otherwise read in count data and normalize
	} else {
		#Get list of samples from directories in rna
		samples <- list.files(path=path1)
		files <- list.files(path=path1,pattern="_ReadsPerGene.out.tab",recursive=TRUE)
		#Get gene names
		genes <- read.table(paste(path1,files[1],sep="/"),header=FALSE,sep="\t",skip=4)[,1]
		#Modify gene names
		genes <- gsub(".Araport11.447","",genes)
		#Read in files and build the table
		df1 <- do.call(cbind,lapply(files,
			function(fn)read.table(paste(path1,fn,sep="/"),header=FALSE,sep="\t",skip=4)[,2]))
		#Add genes as rownames
		row.names(df1)<- genes
		#Add samples as column names 
		colnames(df1) <- samples
		#Create dataframe describing samples & replicates
		coldata <- data.frame(row.names=colnames(df1),condition=gsub("-.$","",samples))
		#Create a DESeqDataSet from the matrix and coldata
		dds <- DESeqDataSetFromMatrix(countData=df1, colData=coldata, design = ~ condition)
		#Normalize data
		dds <- estimateSizeFactors(dds)
		#Extract out normalized counts
		df2 <- data.frame(counts(dds,normalized=TRUE), check.names = F)
		#Create an empty dataframe to get means of different conditions
		df3 <- data.frame(row.names=row.names(df2))
		#Create a new table, where we get the mean of the conditions with replicates
		#Since some conditions only had one replicate, we just keep the original value 
		#for those conditions
		for(x in unique(coldata$condition)){
			if(ncol(data.frame(df2[,grepl(x, colnames(df2))])) > 1){
				df3[,x] <- rowMeans(df2[,grepl(x, colnames(df2))])
			} else {
			df3[,x] <- df2[,grepl(x, colnames(df2))]
			}
		}
	}
	#Output as a table
	write.csv(df3,paste(path2,"/",i,"_normalized_expression.csv",sep=""),
		quote=FALSE,row.names=TRUE)
}
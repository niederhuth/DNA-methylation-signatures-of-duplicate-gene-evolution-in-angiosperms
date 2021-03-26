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

for(i in c("Athaliana")){
	#Input data path
	path1=paste(i,"/rna",sep="")
	#Output data path
	path2=paste("../figures_tables/",i,sep="")
	#Get list of samples from directories in rna
	samples <- list.files(path=path1)
	files <- list.files(path=path1,pattern="_ReadsPerGene.out.tab",recursive=TRUE)
	#Get gene names
	genes <- read.table(paste(path1,files[1],sep="/"),header=FALSE,sep="\t",skip=4)[,1]
	#Modify gene names
	genes <- gsub(".Araport11.447","",genes)
	#Read in files and build the table
	df <- do.call(cbind,lapply(files,
		function(fn)read.table(paste(path1,fn,sep="/"),header=FALSE, sep="\t", skip=4)[,2]))
	#Add genes as rownames
	row.names(df)<- genes
	#Add samples as column names 
	colnames(df) <- samples
	#Output as a table
	write.csv(df,paste(path2,"/",i,"_gene_exp.csv",sep=""),quote=FALSE,row.names=TRUE)
	#Create dataframe describing samples & replicates
	coldata <- data.frame(row.names=colnames(df),condition=gsub("-.$","",samples))
	#Create a DESeqDataSet from the matrix and coldata
	dds <- DESeqDataSetFromMatrix(countData=df, colData=coldata, design = ~ condition)
	#Normalize data
	dds <- estimateSizeFactors(dds)
	#Extract out normalized counts
	df2 <- data.frame(counts(dds,normalized=TRUE))
	#Create an empty dataframe to get means of different conditions
	df3 <- data.frame(row.names=row.names(df2))
	#Create a new table, where we get the mean of the conditions with replicates
	#Since some conditions only had one replicate, we just keep the original value for those conditions
	for(x in unique(coldata$condition)){
		if(ncol(data.frame(df2[,grepl(x, colnames(df2))])) > 1){
			df3[,x] <- rowMeans(df2[,grepl(x, colnames(df2))])
		} else {
			df3[,x] <- df2[,grepl(x, colnames(df2))]
		}
		#We are going to set genes with counts less than 1 (or some other number) to 0
		df3[,x] <- ifelse(df3[,x] < 1,0,df3[,x])
	}
	#We need to add some number to the data, so that we can log transform samples with 0 read counts
	#without getting infinite values
	df3 <- df3 + 1
	#We now create a new data frame, applying the tau tissue specificity 
	df4 <- data.frame(tau=apply(df3,1,tauSpecificity))
	#Add gene names as a column called features
	df4$Feature <- row.names(df4)
	#Read in KaKs table, with duplicate gene pairs and other relevant info
	df5 <- read.csv(paste("../figures_tables/",i,"/",i,"_KaKs_values.csv",sep=""),
		header=TRUE,row.names=1)
	#Merge the data frames together.
	df6 <- merge(df5,df4,by="Feature")
	df7 <- merge(df6,df6,by.x="Duplicate.2",by.y="Feature")
}









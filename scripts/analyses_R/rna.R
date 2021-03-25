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
	coldata <- data.frame(rownames=colnames(df),condition=gsub("-.$","",samples))
	#Create a DESeqDataSet from the matrix and coldata
	dds <- DESeqDataSetFromMatrix(countData=df, colData=coldata, design = ~ condition)
	dds <- estimateSizeFactors(dds)
	dds <- estimateDispersions(dds,fitType="parametric")
}
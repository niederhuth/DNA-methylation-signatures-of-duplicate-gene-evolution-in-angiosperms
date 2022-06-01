library(reshape2)
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

species=c("Sbicolor")

for(i in species){
	path1=paste(i,"/rna/Gene_Expression.tsv")
	df1 <- read.table(path1,header=F,sep="\t")
	df1 <- df1[c(1,3,2)]
	colnames(df1) <- c("Feature","variable","value")
	df2 <- dcast(df1,formula=Feature~variable,fun.aggregate=sum,value.var="value")
	row.names(df2) <- df2$Feature
	df2 <- df2[2:ncol(df2)]
	df3 <- df2
	df3 <- df3 + 1
	df4 <- data.frame(tau=apply(df3,1,tauSpecificity))
	df4$Feature <- row.names(df4)
	#Read in KaKs table, with duplicate gene pairs and other relevant info
	df5 <- read.csv(paste("../figures_tables/",i,"/",i,"_KaKs_values.csv",sep=""),
		header=TRUE,row.names=1)
	#Merge the data frames together.
	df6 <- merge(df5,df4,by="Feature")
} 

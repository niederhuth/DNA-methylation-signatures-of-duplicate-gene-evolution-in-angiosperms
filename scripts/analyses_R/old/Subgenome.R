library(ggplot2)
library(scales) 

species=c("Boleracea","Zmays")

path1 <- "../figures_tables/Subgenome"
if(!file.exists(path1)){
    dir.create(path1)
}

for(a in species){
	path2 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
	path3 <- paste("../misc/",a,"_subgenome.csv",sep="")
	df1 <- read.table(path2,header=TRUE,sep="\t")
	df2 <- read.table(path3,header=TRUE,sep=",")
	df3 <- merge(df2,df1,by="Feature")
	df4 <- data.frame(table(df3$Subgenome,df3$Classification))
	colnames(df4) <- c("Subgenome","Classification","Freq")
	df4 <- rbind(df4,data.frame(Subgenome=c('Total'),
		Classification=data.frame(table(df1$Classification))$Var1,
		Freq=data.frame(table(df1$Classification))$Freq))
	df4$Perc <- NA
	for(b in unique(df4$Subgenome)){
		df4[df4$Subgenome==b,]$Perc <- df4[df4$Subgenome==b,]$Freq/sum(df4[df4$Subgenome==b,]$Freq)
	}

	p <- ggplot(df4[df4$Classification!="Unclassified",]) + 
		geom_bar(aes(x=Classification,y=Perc,fill=Subgenome),
			position="dodge",stat="identity") + 
		theme_bw() + scale_y_continuous("Percentage of genes",
			labels=percent,expand=c(0,0)) +
		theme(axis.text=element_text(color="black"),
			axis.ticks=element_line(color="black"))
	ggsave(paste(path1,"/",a,"_subgenome.pdf",sep=""),p,device="pdf")
}

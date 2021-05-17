library(ggplot2)
library(gghalves)

#Function to calculate tau for expression specificity. 
#See Yanai et al. 2015 Bioinformatics
tauSpecificity <- function(x){
	i = 0
	for(j in c(1:length(x))){
		i = i + 1 - (log2(x[j])/log2(max(x)))
	}
	y = i/(length(x) - 1)
	return(y)
}

species=c("Athaliana","Gmax","Pvulgaris","Sbicolor")

for(i in species){
	#Set output path
	path1 <- paste("../figures_tables/",i,"/",sep="")
	#read in normalized expression table
	df1 <- read.csv(paste(path1,i,"_normalized_expression.csv",sep=""),
		header=TRUE,row.names=1)
	#We will add 1 to the read counts as taking log2 of 0 results in infinity values
	df1 <- df1 + 1
	#We now create a new data frame, applying the tau tissue specificity 
	df2 <- data.frame(tau=apply(df1,1,tauSpecificity))
	#Add gene names as a column called features
	df2$Feature <- row.names(df2)
	#Read in methylation classifications
	df3 <- read.table(paste(i,"/methylpy/results/",i,"_classified_genes.tsv",sep=""),
		sep="\t",header=TRUE)[c("Feature","Classification")]
	colnames(df3) <- c("Feature","Methylation")
	#Read in duplication classifications
	df4 <- read.table(paste(i,"/dupgen/results-unique/classified_genes.tsv",sep=""),
		sep="\t",header=TRUE)
	#Merge df3 & df2
	df5 <- merge(merge(df3,df4,by="Feature"),df2,by="Feature")
	#Write as table
	write.csv(df5,paste(path1,i,"_gene_tau.csv",sep=""),quote=FALSE,row.names=FALSE)
	#Filter and add order
	df6 <- df5[!is.na(df5$Methylation) & df5$Methylation != "Unclassified",]
	df6$order <- ifelse(df6$Methylation=="gbM",1,NA)
	df6$order <- ifelse(df6$Methylation=="teM",2,df6$order)
	df6$order <- ifelse(df6$Methylation=="unM",3,df6$order)
	df6$order2 <- ifelse(df6$Duplication=="wgd",1,NA)
	df6$order2 <- ifelse(df6$Duplication=="tandem",2,df6$order2)
	df6$order2 <- ifelse(df6$Duplication=="proximal",3,df6$order2)
	df6$order2 <- ifelse(df6$Duplication=="translocated",4,df6$order2)
	df6$order2 <- ifelse(df6$Duplication=="dispersed",5,df6$order2)

	#Make Plots of genic methylation & tau
	p <- ggplot(df6) +
		theme_bw() +
		geom_half_boxplot(aes(x=reorder(Methylation,order),y=tau,fill=Methylation),
			side="l",errorbar.draw=TRUE,outlier.color=NA) +
		geom_half_point(aes(x=reorder(Methylation,order),y=tau,color=Methylation),
			side="r",size=0.5) +
		labs(x="",y="Tau")
	ggsave(paste(path1,i,"_Tau_metclass_gghalves.pdf",sep=""),p,device="pdf",width=4,height=3)
	#Make Plots of Duplication & tau
	p <- ggplot(df6[df6$Duplication!="unclassified" & df6$Duplication!="singletons",]) +
		theme_bw() +
		geom_boxplot(aes(x=reorder(Duplication,order2),y=tau,fill=Duplication),
			outlier.color=NA) +
		labs(x="",y="Tau")
	ggsave(paste(path1,i,"_Tau_duplication.pdf",sep=""),p,device="pdf",width=4,height=3)
	#Plot Duplication by Methylation
	#gbM
	p <- ggplot(df6[df6$Duplication!="unclassified" & df6$Duplication!="singletons" &
					df6$Methylation=="gbM",]) +
		theme_bw() +
		geom_boxplot(aes(x=reorder(Duplication,order2),y=tau,fill=Duplication),
			outlier.color=NA) +
		labs(x="",y="Tau")
	ggsave(paste(path1,i,"_Tau_duplication_gbM.pdf",sep=""),p,device="pdf",width=4,height=3)
	#teM
	p <- ggplot(df6[df6$Duplication!="unclassified" & df6$Duplication!="singletons" &
					df6$Methylation=="teM",]) +
		theme_bw() +
		geom_boxplot(aes(x=reorder(Duplication,order2),y=tau,fill=Duplication),
			outlier.color=NA) +
		labs(x="",y="Tau")
	ggsave(paste(path1,i,"_Tau_duplication_teM.pdf",sep=""),p,device="pdf",width=4,height=3)
	#unM
	p <- ggplot(df6[df6$Duplication!="unclassified" & df6$Duplication!="singletons" &
					df6$Methylation=="unM",]) +
		theme_bw() +
		geom_boxplot(aes(x=reorder(Duplication,order2),y=tau,fill=Duplication),
			outlier.color=NA) +
		labs(x="",y="Tau")
	ggsave(paste(path1,i,"_Tau_duplication_unM.pdf",sep=""),p,device="pdf",width=4,height=3)

	#Read in Duplicate pair methylation similarity
	df7 <- read.csv(paste("../figures_tables/",i,"/",i,"_Duplicate_pair_met.csv",sep=""),
		header=TRUE)
	#Merge the dataframes
	df8 <- merge(df7,df2,by.x="Duplicate.1",by.y="Feature")
	df8 <- merge(df8,df2,by.x="Duplicate.2",by.y="Feature")
	df8 <- df8[,c(1,3,2,4:9)]
	colnames(df8) <- c(colnames(df7),"Duplicate.1_tau","Duplicate.2_tau")
	#absolute difference between duplicate pair tau
	df8$absdiff <- abs(df8$Duplicate.1_tau-df8$Duplicate.2_tau)
	#Correlate expression
	df8$Correlation <- NA
	for(a in 1:nrow(df8)){
		df8[a,"Correlation"] <- cor(t(df1[row.names(df1)==df8[a,"Duplicate.1"],]),
									t(df1[row.names(df1)==df8[a,"Duplicate.2"],]))
	}
	#Output df8 as a table
	write.csv(df8,paste(path1,i,"_gene_pair_expression.csv",sep=""),quote=FALSE,row.names=FALSE)
	#Remove unclassified gene pairs for plotting
	df9 <- df8[df8$Duplicate.1_Methylation != "Unclassified" & 
				df8$Duplicate.2_Methylation != "Unclassified",]

	#Plot Correlations
	p <- ggplot(df9) + geom_density(aes(x=Correlation,color=Classification)) +
			theme_bw() + theme()
	ggsave(paste(path1,i,"_expression_correlation_density.pdf",sep=""),p,device="pdf")
	#Absolute difference gghalves
	p <- ggplot(df9) +
			theme_bw() +
			geom_half_boxplot(aes(x=Classification,y=absdiff,fill=Classification),
				side="l",errorbar.draw=TRUE,outlier.color=NA) +
			geom_half_point(aes(x=Classification,y=absdiff,color=Classification), 
				side="r",size=0.5) 
	ggsave(paste(path1,i,"_Tau_AbsDiff_gghalves.pdf",sep=""),p,device="pdf",width=6,height=3)
	#Plot tau of duplicates that switch or do not switch
	#gbM
	df10 <- df9[df9$Duplicate.1_Methylation=="gbM",
				c("Duplicate.1","Duplicate.2_Methylation","Duplicate.1_tau")]
	tmp <- df9[df9$Duplicate.2_Methylation=="gbM",
				c("Duplicate.2","Duplicate.1_Methylation","Duplicate.2_tau")]
	colnames(tmp) <- colnames(df10)
	df10 <- rbind(df10,tmp)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="gbM",2,NA)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="teM",3,df10$order)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="unM",4,df10$order)
	tmp <- df10
	tmp$Duplicate.2_Methylation <- "All gbM"
	tmp$order <- 1
	df10 <- rbind(df10,tmp)
	p <- ggplot(df10) +
			geom_boxplot(aes(reorder(x=Duplicate.2_Methylation,order),
				y=Duplicate.1_tau,fill=Duplicate.2_Methylation)) +
			theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=10), 
				axis.text.y=element_text(color="black",size=10), 
				axis.title.y=element_text(size=12),legend.text=element_text(size=10), 
				panel.border=element_rect(colour="#999999", fill=NA), 
				panel.background=element_rect(fill="white"), 
				panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
				strip.background=element_rect(fill="white"), 
				strip.text.y=element_text(angle=0, hjust=0, size=10),
				plot.title=element_text(size=12),legend.title=element_blank(),
				legend.position="none") + 
			scale_fill_brewer(palette="Paired") +
			labs(x="Methylation of other duplicate",y="Tau - gbM genes  \n")
	ggsave(paste(path1,i,"_Tau_Switching_gbM.pdf",sep=""), p,device="pdf",width=4,height=3)
	#teM
	df10 <- df9[df9$Duplicate.1_Methylation=="teM",
				c("Duplicate.1","Duplicate.2_Methylation","Duplicate.1_tau")]
	tmp <- df9[df9$Duplicate.2_Methylation=="teM",
				c("Duplicate.2","Duplicate.1_Methylation","Duplicate.2_tau")]
	colnames(tmp) <- colnames(df10)
	df10 <- rbind(df10,tmp)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="teM",2,NA)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="gbM",3,df10$order)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="unM",4,df10$order)
	tmp <- df10
	tmp$Duplicate.2_Methylation <- "All teM"
	tmp$order <- 1
	df10 <- rbind(df10,tmp)
	p <- ggplot(df10) +
			geom_boxplot(aes(reorder(x=Duplicate.2_Methylation,order),
				y=Duplicate.1_tau,fill=Duplicate.2_Methylation)) +
			theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=10), 
				axis.text.y=element_text(color="black",size=10), 
				axis.title.y=element_text(size=12),legend.text=element_text(size=10), 
				panel.border=element_rect(colour="#999999", fill=NA), 
				panel.background=element_rect(fill="white"), 
				panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
				strip.background=element_rect(fill="white"), 
				strip.text.y=element_text(angle=0, hjust=0, size=10),
				plot.title=element_text(size=12),legend.title=element_blank(),
				legend.position="none") + 
			scale_fill_brewer(palette="Paired") +
			labs(x="Methylation of other duplicate",y="Tau - teM genes  \n")
	ggsave(paste(path1,i,"_Tau_Switching_teM.pdf",sep=""), p,device="pdf",width=4,height=3)
	#unM
	df10 <- df9[df9$Duplicate.1_Methylation=="unM",
				c("Duplicate.1","Duplicate.2_Methylation","Duplicate.1_tau")]
	tmp <- df9[df9$Duplicate.2_Methylation=="unM",
				c("Duplicate.2","Duplicate.1_Methylation","Duplicate.2_tau")]
	colnames(tmp) <- colnames(df10)
	df10 <- rbind(df10,tmp)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="unM",2,NA)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="gbM",3,df10$order)
	df10$order <- ifelse(df10$Duplicate.2_Methylation=="teM",4,df10$order)
	tmp <- df10
	tmp$Duplicate.2_Methylation <- "All unM"
	tmp$order <- 1
	df10 <- rbind(df10,tmp)
	p <- ggplot(df10) +
			geom_boxplot(aes(reorder(x=Duplicate.2_Methylation,order),
				y=Duplicate.1_tau,fill=Duplicate.2_Methylation)) +
			theme(axis.text.x=element_text(color="black",angle=45,hjust=1,size=10), 
				axis.text.y=element_text(color="black",size=10), 
				axis.title.y=element_text(size=12),legend.text=element_text(size=10), 
				panel.border=element_rect(colour="#999999", fill=NA), 
				panel.background=element_rect(fill="white"), 
				panel.grid.major=element_line(colour="#cccccc", linetype="dotted"), 
				strip.background=element_rect(fill="white"), 
				strip.text.y=element_text(angle=0, hjust=0, size=10),
				plot.title=element_text(size=12),legend.title=element_blank(),
				legend.position="none") + 
			scale_fill_brewer(palette="Paired") +
			labs(x="Methylation of other duplicate",y="Tau - unM genes  \n")
	ggsave(paste(path1,i,"_Tau_Switching_unM.pdf",sep=""), p,device="pdf",width=4,height=3)
}





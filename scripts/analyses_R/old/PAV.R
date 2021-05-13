library(reshape2)
library(ggplot2)
library(scales)


path1 <- "../../figures_tables/PAV"
if(!file.exists(path1)){
    dir.create(path1)
}


FE <- data.frame()
#Zmays
df1 <- read.table("Zmays/methylpy/results/Zmays_classified_genes.tsv",
	header=TRUE,sep="\t")[,c(1,23)]
df2 <- read.csv("../misc/Zmays_PAV.csv",header=TRUE)
df3 <- melt(df2)
colnames(df3) <- c("Feature","Variety","PAV")
df4 <- merge(df3,df1,by="Feature")
df5 <- data.frame(table(unique(df4[df4$PAV<0.2,c(1,4)])$Classification),
	table(df1$Classification))[,c(1,2,4)]
colnames(df5) <- c("Classification","PAV","Genes")
df5 <- rbind(df5,data.frame(Classification=c("Total"),
	PAV=sum(df5$PAV),Genes=nrow(df1)))
df5 <- df5[c(1,2,4,5),]
df5$Percent <- df5$PAV/df5$Genes
df5$order <- c(2,3,4,1)
p <- ggplot(df5) + 
	geom_bar(aes(y=Percent,x=reorder(Classification,order),
		fill=Classification),stat="identity") + 
	scale_y_continuous("Percent PAV Genes",labels=percent,expand=c(0,0)) + 
	theme_bw() + 
	theme(axis.text=element_text(color="black"),
		axis.ticks=element_line(color="black"),
		legend.position="None") + xlab("") + 
	scale_fill_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7"))
ggsave(paste(path1,"/Zmays_PAV.pdf",sep=""),p,device="pdf",height=4,width=4)
for(i in 1:3){
	FE <- rbind(FE,data.frame(species=c("Z. mays"),classification=df5[i,1],
	estimate=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$estimate,
	p.value=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$p.value))
}

#Boleracea: https://appliedbioinformatics.com.au/download/BOLEPan.pav.13062016.vcf.gz
df1 <- read.table("Boleracea/methylpy/results/Boleracea_classified_genes.tsv",
	header=TRUE,sep="\t")[,c(1,23)]
df2 <- read.table("../../data/BOLEPan.pav.13062016.vcf",header=FALSE,
	sep="\t")[,c(3,10:19)]
colnames(df2) <- c("Feature","Broccoli","Brussels","Cabbage1","Cabbage2",
	"Cauliflower1","Cauliflower2","Kale","Kohohlrabi","TO1000","Macrocarpa")
for(i in colnames(df2[2:length(df2)])){
	df2[,i] <- gsub("1/1",1,df2[,i])
	df2[,i] <- gsub("0/0",0,df2[,i])
	df2[,i] <- as.numeric(df2[,i])
}
df3 <- melt(df2)
colnames(df3) <- c("Feature","Variety","PAV")
df4 <- merge(df3,df1,by="Feature")
df5 <- data.frame(table(unique(df4[df4$PAV==0,c(1,4)])$Classification),
	table(df1$Classification))[,c(1,2,4)]
colnames(df5) <- c("Classification","PAV","Genes")
df5 <- rbind(df5,data.frame(Classification=c("Total"),
	PAV=sum(df5$PAV),Genes=nrow(df1)))
df5 <- df5[c(1,2,4,5),]
df5$Percent <- df5$PAV/df5$Genes
df5$order <- c(2,3,4,1)
p <- ggplot(df5) + 
	geom_bar(aes(y=Percent,x=reorder(Classification,order),
		fill=Classification),stat="identity") + 
	scale_y_continuous("Percent PAV Genes",labels=percent,expand=c(0,0)) + 
	theme_bw() + 
	theme(axis.text=element_text(color="black"),
		axis.ticks=element_line(color="black"),
		legend.position="None") + xlab("") + 
	scale_fill_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7"))
ggsave(paste(path1,"/Boleracea_PAV.pdf",sep=""),p,device="pdf",height=4,width=4)
for(i in 1:3){
	FE <- rbind(FE,data.frame(species=c("B. oleracea"),classification=df5[i,1],
	estimate=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$estimate,
	p.value=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$p.value))
}

#Slycopersicum #https://www.nature.com/articles/s41588-019-0410-2#Sec23
df1 <- read.table("Slycopersicum/methylpy/results/Slycopersicum_classified_genes.tsv",
	header=TRUE,sep="\t")[,c(1,23)]
df2 <- read.table("../../data/Slycopersicum_PAV.txt",header=TRUE,
	sep="\t")
df4 <- merge(df2,df1,by="Feature")
df5 <- data.frame(table(df4$Classification),
	table(df1$Classification))[,c(1,2,4)]
colnames(df5) <- c("Classification","PAV","Genes")
df5 <- rbind(df5,data.frame(Classification=c("Total"),
	PAV=sum(df5$PAV),Genes=nrow(df1)))
df5 <- df5[c(1,2,4,5),]
df5$Percent <- df5$PAV/df5$Genes
df5$order <- c(2,3,4,1)
p <- ggplot(df5) + 
	geom_bar(aes(y=Percent,x=reorder(Classification,order),
		fill=Classification),stat="identity") + 
	scale_y_continuous("Percent PAV Genes",labels=percent,expand=c(0,0)) + 
	theme_bw() + 
	theme(axis.text=element_text(color="black"),
		axis.ticks=element_line(color="black"),
		legend.position="None") + xlab("") + 
	scale_fill_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7"))
ggsave(paste(path1,"/Slycopersicum_PAV.pdf",sep=""),p,device="pdf",height=4,width=4)
for(i in 1:3){
	FE <- rbind(FE,data.frame(species=c("S. lycoperiscum"),classification=df5[i,1],
	estimate=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$estimate,
	p.value=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$p.value))
}

#Stuberosum http://www.plantcell.org/highwire/filestream/4514/field_highwire_adjunct_files/7/TPC2015-00538-RAR2_Supplemental_Data_Set_5.xlsx
df1 <- read.table("Stuberosum/methylpy/results/Stuberosum_classified_genes.tsv",
	header=TRUE,sep="\t")[,c(1,23)]
df2 <- read.csv("../../data/TPC2015-00538-RAR2_Supplemental_Data_Set_5.csv",
	header=TRUE)[,c(1,6:17)]
df3 <- melt(df2)
colnames(df3) <- c("Feature","Variety","PAV")
df4 <- merge(df3,df1,by="Feature")
df5 <- data.frame(table(unique(df4[df4$PAV<0.2,c(1,4)])$Classification),
	table(df1$Classification))[,c(1,2,4)]
colnames(df5) <- c("Classification","PAV","Genes")
df5 <- rbind(df5,data.frame(Classification=c("Total"),
	PAV=sum(df5$PAV),Genes=nrow(df1)))
df5 <- df5[c(1,2,4,5),]
df5$Percent <- df5$PAV/df5$Genes
df5$order <- c(2,3,4,1)
p <- ggplot(df5) + 
	geom_bar(aes(y=Percent,x=reorder(Classification,order),
		fill=Classification),stat="identity") + 
	scale_y_continuous("Percent PAV Genes",labels=percent,expand=c(0,0)) + 
	theme_bw() + 
	theme(axis.text=element_text(color="black"),
		axis.ticks=element_line(color="black"),
		legend.position="None") + xlab("") + 
	scale_fill_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7"))
ggsave(paste(path1,"/Stuberosum_PAV.pdf",sep=""),p,device="pdf",height=4,width=4)
for(i in 1:3){
	FE <- rbind(FE,data.frame(species=c("S. tuberosum"),classification=df5[i,1],
	estimate=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$estimate,
	p.value=fisher.test(matrix(c(df5[i,2],
		df5[i,3]-df5[i,2],
		df5[4,2]-df5[i,2],
		df5[4,3]-df5[4,2]-df5[i,3]),
		nrow=2,ncol=2),alternative="two.sided")$p.value))
}

FE$p.adjust <- p.adjust(FE$p.value,method="BH")
write.csv(FE,paste(path1,"/PAV_Fishers_Exact.csv",sep=""),quote=FALSE,row.names=FALSE)
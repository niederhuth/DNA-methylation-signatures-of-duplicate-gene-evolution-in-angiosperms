library(reshape2)
library(ggplot2)
library(scales)

#List species to be analyzed
species <- data.frame(Species=c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
	"Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
	"Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
	"Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
	"Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
	"Vvinifera","Zmays"),Order=c(35,34,20,19,1,5,23,22,40,14,17,26,27,18,28,13,2,21,32,33,37,
	16,38,3,31,24,43,39,4,11,29,25,10,36,30,7,8,41,42,9,15,12,6))

#setwd("data/")
df12 <- df11 <- data.frame()

for(a in species$Species){
	path1 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
	path2 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
	path3 <- paste("../figures_tables/",a,"/",sep="")
	if(!file.exists(path3)){
		dir.create(path3)
	}
	df1 <- read.table(path1,header=TRUE,sep="\t")
	df2 <- read.table(path2,header=TRUE,sep="\t")
	df2$Classification <- gsub("Unmethylated","unM",df2$Classification)

	df7 <- data.frame()

	for(b in c("wgd","proximal","dispersed","tandem","transposed")){
		path4u <- paste(a,"/dupgen/results-unique/",a,".",b,".pairs-unique",sep="")
		path4a <- paste(a,"/dupgen/results/",a,".",b,".pairs",sep="")
		df3u <- read.table(path4u,header=TRUE,sep="\t")[,c(1,3,5)]
		df3a <- read.table(path4a,header=TRUE,sep="\t")[,c(1,3,5)]
		colnames(df3u) <- colnames(df3a) <- c("Duplicate.1","Duplicate.2","E.value")
		df4 <- data.frame()
		for(gene in df1[df1$Duplication==b,]$Feature){
			tmp <- df3u[df3u$Duplicate.1==gene | df3u$Duplicate.2==gene,]
			if(nrow(tmp) == 0){
				tmp <- df3a[df3a$Duplicate.1==gene | df3a$Duplicate.2==gene,]
			}
			df4 <- rbind(df4,tmp[row.names(tmp) == sample(row.names(tmp[tmp$E.value==min(tmp$E.value),]),1),c(1,2)])

		}
		df5 <- merge(df4,df2[,c(1,30)],by.x="Duplicate.1",by.y="Feature")
		df6 <- na.omit(merge(df5,df2[,c(1,30)],by.x="Duplicate.2",by.y="Feature")[,c(2,3,1,4)])
		colnames(df6) <- c("Duplicate.1","Duplicate.1_Methylation","Duplicate.2","Duplicate.2_Methylation")
		df6$Duplication <- c(b)
		df6$Similarity <- ifelse(df6$Duplicate.1_Methylation==df6$Duplicate.2_Methylation,"Same","Different")
		df6$Classification <- "NA"

		for(row in 1:nrow(df6)){
			if(df6[row,]$Duplicate.1_Methylation == df6[row,]$Duplicate.2_Methylation){
				df6[row,]$Classification = paste(df6[row,]$Duplicate.1_Methylation,df6[row,]$Duplicate.2_Methylation,sep="-")
			} else if(df6[row,]$Duplicate.1_Methylation == "gbM" & df6[row,]$Duplicate.2_Methylation == "Unclassified" || df6[row,]$Duplicate.1_Methylation == "Unclassified" & df6[row,]$Duplicate.2_Methylation == "gbM" ){
				df6[row,]$Classification = "gbM-Unclassified"
			}  else if(df6[row,]$Duplicate.1_Methylation == "gbM" & df6[row,]$Duplicate.2_Methylation == "teM" || df6[row,]$Duplicate.1_Methylation == "teM" & df6[row,]$Duplicate.2_Methylation == "gbM" ){
				df6[row,]$Classification = "gbM-teM"
			}  else if(df6[row,]$Duplicate.1_Methylation == "gbM" & df6[row,]$Duplicate.2_Methylation == "unM" || df6[row,]$Duplicate.1_Methylation == "unM" & df6[row,]$Duplicate.2_Methylation == "gbM" ){
				df6[row,]$Classification = "gbM-unM"
			} else if(df6[row,]$Duplicate.1_Methylation == "unM" & df6[row,]$Duplicate.2_Methylation == "teM" || df6[row,]$Duplicate.1_Methylation == "teM" & df6[row,]$Duplicate.2_Methylation == "unM" ){
				df6[row,]$Classification = "unM-teM"
			} else if(df6[row,]$Duplicate.1_Methylation == "unM" & df6[row,]$Duplicate.2_Methylation == "Unclassified" || df6[row,]$Duplicate.1_Methylation == "Unclassified" & df6[row,]$Duplicate.2_Methylation == "unM" ){
				df6[row,]$Classification = "unM-Unclassified"
			} else if(df6[row,]$Duplicate.1_Methylation == "teM" & df6[row,]$Duplicate.2_Methylation == "Unclassified" || df6[row,]$Duplicate.1_Methylation == "Unclassified" & df6[row,]$Duplicate.2_Methylation == "teM" ){
				df6[row,]$Classification = "teM-Unclassified"
			} else {
				df6[row,]$Classification = df6[row,]$Classification
			}
		}

		df7 <- rbind(df7,df6)
	}

	df7 <- unique(df7)
	df7$Similarity <- ifelse(df7$Duplicate.1_Methylation == "Unclassified" | df7$Duplicate.2_Methylation == "Unclassified",
		"Undetermined",df7$Similarity)
	df8 <- df7[df7$Duplicate.1_Methylation == df7$Duplicate.2_Methylation,]
	for(x in unique(df7[df7$Duplicate.1_Methylation != df7$Duplicate.2_Methylation,]$Classification)){
		y <- gsub("-.*","",x)
		df8 <- rbind(df8,df7[df7$Classification == x & df7$Duplicate.1_Methylation == y,])
		tmp <- df7[df7$Classification == x & df7$Duplicate.2_Methylation == y,c(3,4,1,2,5,6,7)]
		colnames(tmp) <- colnames(df8)
		df8 <- rbind(df8,tmp)
	}
	df8 <- unique(df8)

	df9 <- data.frame(species=a,table(df8$Duplication,df8$Similarity),PhyloOrder=species[species$Species==a,]$Order)
	colnames(df9) <- c("Species","Duplication","Similarity","Number","PhyloOrder")
	df10 <- data.frame(species=a,table(df8$Duplication,df8$Classification),PhyloOrder=species[species$Species==a,]$Order)
	colnames(df10) <- c("Species","Duplication","Classification","Number","PhyloOrder")
	tmp <-  as.vector(table(df8$Classification))
	names(tmp) <- names(table(df7$Classification))

	write.csv(df8,paste(path3,a,"_duplicate_pair_met.csv",sep=""),quote=FALSE,row.names=FALSE)
	write.csv(df9,paste(path3,a,"_duplicate_similarity_1.csv",sep=""),row.names=FALSE,quote=FALSE)
	write.csv(df10,paste(path3,a,"_duplicate_similarity_2.csv",sep=""),row.names=FALSE,quote=FALSE)
	df11 <- rbind(df11,df9)
	df12 <- rbind(df12,data.frame(species=a,t(tmp)))
}

#Output combined tables
path5 <- paste("../figures_tables/duplicate_similarity")
if(!file.exists(path5)){
	dir.create(path1)
}
write.csv(df11,paste(path5,"/All_DuplicateSimilarity.csv",sep=""),row.names=FALSE,quote=FALSE)
colnames(df12) <- gsub("\\.","-",colnames(df12))
write.csv(df12,paste(path5,"/All_DuplicatePair_classification.csv",sep=""),row.names=FALSE,quote=FALSE)

#### Duplicate-Similarity figures ####
df11$Percent <- NA
for(i in 1:nrow(df9)){
	df11[i,]$Percent <- df11[i,]$Number/
	sum(df11[df9$Species==df11[i,]$Species & df11$Duplication==df11[i,]$Duplication,]$Number)
}

# To get one figure for all different types of duplicates
p <- ggplot(df11, aes(fill=Similarity,y=Number,x=reorder(Species,PhyloOrder))) + 
	geom_bar(position="stack",stat="identity") +   #position=fill gives percentage, position=stack gives stacked number barplot
	theme_bw() +
	scale_y_continuous("Percentage of genes",expand=c(0,0),labels=percent) + #y-axis shows % rather than numbers
	scale_x_discrete("")
p + facet_wrap(~Duplication, scales="free", nrow =5)
ggsave(paste(path5,"/duplicate_similarity.pdf",sep=""),p,device="pdf")


# To plot each types of duplicates seperately
for(a in c("wgd","proximal","dispersed","tandem","transposed")){
	p <- ggplot(df11[df11$Duplication==a,], aes(fill=Similarity,y=Percent,x=reorder(Species,PhyloOrder))) + 
		geom_bar(position="fill",stat="identity") +   #position=fill gives percentage, position=stack gives stacked number barplot
		theme_bw() +
		scale_y_continuous("Percentage of genes",expand=c(0,0),labels=percent) #y-axis shows % rather than numbers
	ggsave(paste(path5,"/",a,"_duplicate_simiarlity.pdf",sep=""),p,device="pdf")
}

 


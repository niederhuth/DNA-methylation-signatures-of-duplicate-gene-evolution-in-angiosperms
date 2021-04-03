library(reshape2)
library(ggplot2)

#List species to be analyzed
species = c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
	"Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
	"Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
	"Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
	"Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
	"Vvinifera","Zmays")

setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/data/")

for(a in species){
	path1 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
	path2 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
	path3 <- paste("../figures_tables/",a,"/",sep="")
	df1 <- read.table(path1,header=TRUE,sep="\t")
	df2 <- read.table(path2,header=TRUE,sep="\t")
	df2$Classification <- gsub("Unmethylated","unM",df2$Classification)

	df6 <- data.frame()

	for(b in c("wgd","proximal","dispersed","tandem","transposed")){
		path4 <- paste(a,"/dupgen/results-unique/",a,".",b,".pairs-unique",sep="")
		df3 <- read.table(path4,header=TRUE,sep="\t")[,c(1,3)]
		colnames(df3) <- c("Duplicate.1","Duplicate.2")
		df4 <- merge(df3,df2[,c(1,30)],by.x="Duplicate.1",by.y="Feature")
		df5 <- na.omit(merge(df4,df2[,c(1,30)],by.x="Duplicate.2",by.y="Feature")[,c(2,3,1,4)])
		colnames(df5) <- c("Duplicate.1","Duplicate.1_Methylation","Duplicate.2","Duplicate.2_Methylation")
		df5$Duplication <- c(b)
		df5$Similarity <- ifelse(df5$Duplicate.1_Methylation==df5$Duplicate.2_Methylation,"Same","Different")
		df5$Classification <- "NA"

		for(row in 1:nrow(df5)){
			if(df5[row,]$Duplicate.1_Methylation == df5[row,]$Duplicate.2_Methylation){
				df5[row,]$Classification = paste(df5[row,]$Duplicate.1_Methylation,df5[row,]$Duplicate.2_Methylation,sep="-")
			} else if(df5[row,]$Duplicate.1_Methylation == "gbM" & df5[row,]$Duplicate.2_Methylation == "Unclassified" || df5[row,]$Duplicate.1_Methylation == "Unclassified" & df5[row,]$Duplicate.2_Methylation == "gbM" ){
				df5[row,]$Classification = "gbM-Unclassified"
			}  else if(df5[row,]$Duplicate.1_Methylation == "gbM" & df5[row,]$Duplicate.2_Methylation == "teM" || df5[row,]$Duplicate.1_Methylation == "teM" & df5[row,]$Duplicate.2_Methylation == "gbM" ){
				df5[row,]$Classification = "gbM-teM"
			}  else if(df5[row,]$Duplicate.1_Methylation == "gbM" & df5[row,]$Duplicate.2_Methylation == "unM" || df5[row,]$Duplicate.1_Methylation == "unM" & df5[row,]$Duplicate.2_Methylation == "gbM" ){
				df5[row,]$Classification = "gbM-unM"
			} else if(df5[row,]$Duplicate.1_Methylation == "unM" & df5[row,]$Duplicate.2_Methylation == "teM" || df5[row,]$Duplicate.1_Methylation == "teM" & df5[row,]$Duplicate.2_Methylation == "unM" ){
				df5[row,]$Classification = "unM-teM"
			} else if(df5[row,]$Duplicate.1_Methylation == "unM" & df5[row,]$Duplicate.2_Methylation == "Unclassified" || df5[row,]$Duplicate.1_Methylation == "Unclassified" & df5[row,]$Duplicate.2_Methylation == "unM" ){
				df5[row,]$Classification = "unM-Unclassified"
			} else if(df5[row,]$Duplicate.1_Methylation == "teM" & df5[row,]$Duplicate.2_Methylation == "Unclassified" || df5[row,]$Duplicate.1_Methylation == "Unclassified" & df5[row,]$Duplicate.2_Methylation == "teM" ){
				df5[row,]$Classification = "teM-Unclassified"
			} else {
				df5[row,]$Classification = df5[row,]$Classification
			}
		}

		df6 <- rbind(df6,df5)
	}


	df6$Similarity <- ifelse(df6$Duplicate.1_Methylation == "Unclassified" | df6$Duplicate.2_Methylation == "Unclassified",
		"Undetermined",df6$Similarity)
	df7 <- df6[df6$Duplicate.1_Methylation == df6$Duplicate.2_Methylation,]
	for(x in unique(df6[df6$Duplicate.1_Methylation != df6$Duplicate.2_Methylation,]$Classification)){
		y <- gsub("-.*","",x)
		df7 <- rbind(df7,df6[df6$Classification == x & df6$Duplicate.1_Methylation == y,])
		tmp <- df6[df6$Classification == x & df6$Duplicate.2_Methylation == y,c(3,4,1,2,5,6,7)]
		colnames(tmp) <- colnames(df7)
		df7 <- rbind(df7,tmp)
	}
	df7 <- unique(df7)

	df8 <- data.frame(species=a,table(df7$Duplication,df7$Similarity))
	colnames(df8) <- c("Species","Duplication","Similarity","Number")
	df9 <- data.frame(species=a,table(df7$Duplication,df7$Classification))
	colnames(df9) <- c("Species","Duplication","Classification","Number")

	write.csv(df7,paste(path3,a,"_Duplicate_pair_met.csv",sep=""),quote=FALSE,row.names=FALSE)
	write.csv(df8,paste(path3,a,"_duplicate_similarity_1.csv",sep=""),row.names=FALSE,quote=FALSE)
	write.csv(df9,paste(path3,a,"_duplicate_similarity_2.csv",sep=""),row.names=FALSE,quote=FALSE)
}


#### Duplicate-Similarity figure ####
​
# First to make a file for all species 
setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/figures_tables/")
​
df9 <- read.csv("All_DuplicateSimilarity.csv") # this could be improved, but for now, I create an empty 
# ouptut file and add all files into this
for (a in species){
  path <- paste(a,"/", a,"_Duplicate_Similarity_1.csv",sep="")
  df10 <- read.csv(path, header=TRUE)
  df9 <- rbind(df9, df10)
}
write.csv(df9, file ="All_DuplicateSimilarity.csv")
​
# Using a subset datset to plot one graph
​
data2 <- read.csv("Dispersed_DuplicateSimilarity.csv")
data2$Similarity <- factor(data2$Similarity, levels(data2$Similarity)[c(3,1,2)])
ggplot(data2, aes(fill=Similarity, y=Number, x=Species)) + 
  geom_bar(position="stack", stat="identity")
​
# Using the whole datset to generate faceted graphs
​
library(scales) # necesssary for percent
​
df11 <- read.csv("All_DuplicateSimilarity.csv")
df11 <- merge(df11,unique(df11[c(1,2)]),by="Species") # for ordering based on slno.
​
df11$Similarity <- factor(df11$Similarity, levels(df11$Similarity)[c(3,1,2)])
​
# To get one figure for all different types of duplicates
p <- ggplot(df12, aes(fill=Similarity, y=Number, x=reorder(Species, slno.x))) + 
  geom_bar(position="fill", stat="identity") +   #position=fill gives percentage, position=stack gives stacked number barplot
  theme_bw() +
  scale_y_continuous("Percentage of genes",expand=c(0,0),labels=percent) #y-axis shows % rather than numbers
p + facet_wrap(~Duplication, scales="free", nrow =5)
​
# To plot each types of duplicates seperately
​
df12 <- subset(df11, Duplication=="wgd")
 q <- ggplot(df12, aes(fill=Similarity, y=Number, x=reorder(Species, slno.x))) + 
    geom_bar(position="fill", stat="identity") +   #position=fill gives percentage, position=stack gives stacked number barplot
    theme_bw() +
   scale_y_continuous("Percentage of genes",expand=c(0,0),labels=percent) #y-axis shows % rather than numbers
 q
 


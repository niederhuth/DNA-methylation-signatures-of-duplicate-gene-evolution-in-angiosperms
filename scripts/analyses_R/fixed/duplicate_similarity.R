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

	df9 <- data.frame(species=a,table(df8$Duplication,df8$Similarity))
	colnames(df9) <- c("Species","Duplication","Similarity","Number")
	df10 <- data.frame(species=a,table(df8$Duplication,df8$Classification))
	colnames(df10) <- c("Species","Duplication","Classification","Number")

	write.csv(df8,paste(path3,a,"_Duplicate_pair_met2.csv",sep=""),quote=FALSE,row.names=FALSE)
	write.csv(df9,paste(path3,a,"_duplicate_similarity2_1.csv",sep=""),row.names=FALSE,quote=FALSE)
	write.csv(df10,paste(path3,a,"_duplicate_similarity2_2.csv",sep=""),row.names=FALSE,quote=FALSE)
}


#### Duplicate-Similarity figure ####

# First to make a file for all species 
setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/figures_tables/")

df10 <- read.csv("All_DuplicateSimilarity2.csv") # this could be improved, but for now, I create an empty 
# ouptut file and add all files into this
for (a in species){
  path <- paste(a,"/", a,"_Duplicate_Similarity2_1.csv",sep="")
  df11 <- read.csv(path, header=TRUE)
  df10 <- rbind(df10, df11)
}
write.csv(df10, file ="All_DuplicateSimilarity2.csv")


# Using the whole datset to generate faceted graphs

library(scales) # necesssary for percent

#All_DuplicateSimilarity2.csv is manually changed to include species order.

df12 <- read.csv("All_DuplicateSimilarity2.csv")
df12 <- merge(df12,unique(df12[c(1,2)]),by="Species") # for ordering based on slno.

df12$Similarity <- factor(df12$Similarity, levels(df12$Similarity)[c(3,1,2)])

# To get one figure for all different types of duplicates
p <- ggplot(df12, aes(fill=Similarity, y=Number, x=reorder(Species, Slno.x))) + 
  geom_bar(position="stack", stat="identity") +   #position=fill gives percentage, position=stack gives stacked number barplot
  theme_bw() +
  scale_y_continuous("Percentage of genes",expand=c(0,0),labels=percent) #y-axis shows % rather than numbers
p + facet_wrap(~Duplication, scales="free", nrow =5)

# To plot each types of duplicates seperately

df13 <- subset(df12, Duplication=="transposed") #Change duplication time for each plot
 q <- ggplot(df13, aes(fill=Similarity, y=Number, x=reorder(Species, Slno.x))) + 
    geom_bar(position="fill", stat="identity") +   #position=fill gives percentage, position=stack gives stacked number barplot
    theme_bw() +
   scale_y_continuous("Percentage of genes",expand=c(0,0),labels=percent) #y-axis shows % rather than numbers
 q
 


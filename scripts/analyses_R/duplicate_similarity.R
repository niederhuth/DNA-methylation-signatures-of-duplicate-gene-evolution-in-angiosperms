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



df8 <- df7 <- data.frame()

for( a in species){
	path1 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
	path2 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
	df1 <- read.table(path1,header=TRUE,sep="\t")
	df2 <- read.table(path2,header=TRUE,sep="\t")

	df6 <- data.frame()

	for(b in c("wgd","proximal","dispersed","tandem","transposed")){
		path3 <- paste(a,"/dupgen/results-unique/",a,".",b,".pairs-unique",sep="")
		df3 <- read.table(path3,header=TRUE,sep="\t")[,c(1,3)]
		colnames(df3) <- c("Duplicate.1","Duplicate.2")
		df4 <- merge(df3,df2[,c(1,23)],by.x="Duplicate.1",by.y="Feature")
		df5 <- na.omit(merge(df4,df2[,c(1,23)],by.x="Duplicate.2",by.y="Feature")[,c(2,3,1,4)])
		colnames(df5) <- c("Duplicate.1","Duplicate.1_Methylation","Duplicate.2","Duplicate.2_Methylation")
		df5$Duplication <- c(b)
		df5$Similarity <- ifelse(df5$Duplicate.1_Methylation==df5$Duplicate.2_Methylation,"Same","Different")
		df5$Classification <- "NA"

		for(row in 1:nrow(df5)){
			if(df5[row,]$Duplicate.1_Methylation == df5[row,]$Duplicate.2_Methylation){
				df5[row,]$Classification = paste(df5[row,]$Duplicate.1_Methylation,df5[row,]$Duplicate.2_Methylation,sep="-")
			} else if(df5[row,]$Duplicate.1_Methylation == "gbM" & df5[row,]$Duplicate.2_Methylation == "Unclassified" || df5[row,]$Duplicate.1_Methylation == "Unclassified" & df5[row,]$Duplicate.2_Methylation == "gbM" ){
				df5[row,]$Classification = "gbM-Unclassified"
			}  else if(df5[row,]$Duplicate.1_Methylation == "gbM" & df5[row,]$Duplicate.2_Methylation == "TE-like" || df5[row,]$Duplicate.1_Methylation == "TE-like" & df5[row,]$Duplicate.2_Methylation == "gbM" ){
				df5[row,]$Classification = "gbM-TE-like"
			}  else if(df5[row,]$Duplicate.1_Methylation == "gbM" & df5[row,]$Duplicate.2_Methylation == "Unmethylated" || df5[row,]$Duplicate.1_Methylation == "Unmethylated" & df5[row,]$Duplicate.2_Methylation == "gbM" ){
				df5[row,]$Classification = "gbM-Unmethylated"
			} else if(df5[row,]$Duplicate.1_Methylation == "Unmethylated" & df5[row,]$Duplicate.2_Methylation == "TE-like" || df5[row,]$Duplicate.1_Methylation == "TE-like" & df5[row,]$Duplicate.2_Methylation == "Unmethylated" ){
				df5[row,]$Classification = "Unmethylated-TE-like"
			} else if(df5[row,]$Duplicate.1_Methylation == "Unmethylated" & df5[row,]$Duplicate.2_Methylation == "Unclassified" || df5[row,]$Duplicate.1_Methylation == "Unclassified" & df5[row,]$Duplicate.2_Methylation == "Unmethylated" ){
				df5[row,]$Classification = "Unmethylated-Unclassified"
			} else if(df5[row,]$Duplicate.1_Methylation == "TE-like" & df5[row,]$Duplicate.2_Methylation == "Unclassified" || df5[row,]$Duplicate.1_Methylation == "Unclassified" & df5[row,]$Duplicate.2_Methylation == "TE-like" ){
				df5[row,]$Classification = "TE-like-Unclassified"
			} else {
				df5[row,]$Classification = df5[row,]$Classification
			}
		}

		df6 <- rbind(df6,df5)
	}

	df6$Similarity <- ifelse(df6$Duplicate.1_Methylation == "Unclassified" | df6$Duplicate.2_Methylation == "Unclassified",
		"Undetermined",df6$Similarity)

	path4 <- paste(a,"/dupgen/results/kaks_results/",a,sep="")

	wgd <- read.table(paste(path4,".wgd.kaks",sep=""),header=T,sep="\t")
	wgd$Duplication <- c("wgd")
	tandem <- read.table(paste(path4,".tandem.kaks",sep=""),header=T,sep="\t")
	tandem$Duplication <- c("tandem")
	proximal <- read.table(paste(path4,".proximal.kaks",sep=""),header=T,sep="\t")
	proximal$Duplication <- c("proximal")
	transposed <- read.table(paste(path4,".transposed.kaks",sep=""),header=T,sep="\t")
	transposed$Duplication <- c("transposed")
	dispersed <- read.table(paste(path4,".dispersed.kaks",sep=""),header=T,sep="\t")
	dispersed$Duplication <- c("dispersed")

	df7 <- rbind(wgd,tandem,proximal,transposed,dispersed)
	df8 <- merge(df6,df7[,c(1:6)],by.x=c('Duplicate.1','Duplicate.2'),
		by.y=c('Duplicate.1','Duplicate.2'))

	df9 <- data.frame(species=a,table(df6$Duplication,df6$Similarity))
	colnames(df9) <- c("Species","Duplication","Similarity","Number")
	df10 <- data.frame(species=a,table(df6$Duplication,df6$Classification))
	colnames(df10) <- c("Species","Duplication","Classification","Number")
}

write.csv(df9,"../figures_tables/duplicate_similarity_1.csv",quote=FALSE)
write.csv(df10,"../figures_tables/duplicate_similarity_2.csv",quote=FALSE)





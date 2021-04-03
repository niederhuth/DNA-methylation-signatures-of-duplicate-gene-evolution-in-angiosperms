
df9 <- df8 <- data.frame()
for(i in c("Athaliana")){
	path1 <- paste(i,"/methylpy/results/",i,sep="")
	path2 <- paste(i,"/dupgen/results-unique/",sep="")
	path3 <- paste("../figures_tables/",i,"/",sep="")
	df1 <- read.table(paste(path1,"_classified_genes.tsv",sep=""),
		header=TRUE,sep="\t")[c("Feature","Classification")]
	df1 <- df1[df1$Classification != "Unclassified",]
	df2 <- read.table(paste(path2,"classified_genes.tsv",sep=""),header=TRUE,sep="\t")
	df2 <- df2[df2$Duplication != "unclassified" & df2$Duplication != "singletons",]
	df3 <- read.table(paste(path1,"_TE_intersections.tsv",sep=""),header=TRUE,sep="\t")
	df3[is.na(df3)] <- 0
	#
	df4 <- merge(df1,df2,by="Feature")
	df5 <- merge(df4,df3[c("Gene","gene_body","up_3000bp","down_3000bp")],
		by.x="Feature",by.y="Gene")
	df5$Classification <- gsub("Unmethylated","unM",df5$Classification)
	df5$Classification <- factor(df5$Classification,levels=c("gbM","unM","teM"))
	df5$Duplication <- factor(df5$Duplication,
		levels=c("wgd","tandem","proximal","transposed","dispersed"))
	df5$TE <- rowSums(df5[c("gene_body","up_3000bp","down_3000bp")])
	df5$TE2 <- ifelse(df5$TE > 0,1,0)
	df6 <- read.csv(paste(path3,i,"_Duplicate_pair_met.csv",sep=""),header=TRUE)
	df6 <- df6[df6$Similarity != "Undetermined",]
	df7 <- merge(df6,df5[c(1,4:8)],by.x="Duplicate.1",by.y="Feature",all.x=TRUE)
	df7 <- merge(df7,df5[c(1,4:8)],by.x="Duplicate.2",by.y="Feature",all.x=TRUE)[c(2,3,1,4:17)]
	df8 <- rbind(df8,data.frame(Species=i,data.frame(table(df5$Classification,df5$TE2))))
	for(x in unique(df7$Classification)){
		df9 <- rbind(df9,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
			Duplicate.2=gsub(".*-","",x),table(df7[df7$Classification == x,c(12,17)])))
	}
}
colnames(df8) <- c("Species","Classification","TE_Presence","Gene_Count")
colnames(df9) <- c("Species","Classification","Duplicate.1","Duplicate.2"
	,"Duplicate.1_TE_presence","Duplicate.2_TE_presence","Gene_Count")
df10 <- df9
#df10 will have all the 1-0 genes for gbM-gbM, unM-unM, and teM-teM set to 0-1
df10[df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==0 & 
	df10$Duplicate.2_TE_presence==1,]$Gene_Count <- df10[df10$Duplicate.1==df10$Duplicate.2 & 
		df10$Duplicate.1_TE_presence==0 & df10$Duplicate.2_TE_presence==1,]$Gene_Count + 
			df10[df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==1 & 
				df10$Duplicate.2_TE_presence==0,]$Gene_Count
#df10 <- df10[!(df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==1 & 
#					df10$Duplicate.2_TE_presence==0),]
df10[df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==1 & 
	df10$Duplicate.2_TE_presence==0,]$Gene_Count <- 0
df10$Gene_Percent <- 0
for(x in 1:nrow(df10)){
	df10[x,]$Gene_Percent <- df10[x,]$Gene_Count/sum(df10[df10$Species==df10[x,]$Species &
																df10$Classification==df10[x,]$Classification,]$Gene_Count)
}
df10$TE_Group <- paste(df10$Duplicate.1_TE_presence,df10$Duplicate.2_TE_presence,sep="-")

ggplot(df10) + geom_bar(aes(x=Classification,y=Gene_Percent,fill=TE_Group),stat="identity")



#Sunil plot code
library(ggplot2)
library(gghalves)
ggplot() +
  theme_bw() +
  scale_y_continuous(expand=c(0,0), limits = c(0,20)) + #points go beyond 10 but Im restricting for better visualization
  geom_half_boxplot(
    data = df4, 
    aes(x = Classification, y = down_10000bp, fill = Classification), 
    side = "l", errorbar.draw = TRUE, # I keep changing the 'y' each time
    outlier.color = NA) +
  geom_half_point(
    data = df4, 
    aes(x = Classification, y = down_10000bp, fill = Classification, color = Classification), 
    side = "r", size =0.5) +
  facet_wrap(~Duplication)
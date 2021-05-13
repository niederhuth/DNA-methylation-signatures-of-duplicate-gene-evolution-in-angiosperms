
df17 <- df11 <- df9 <- df8 <- data.frame()
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
	colnames(df8) <- c("Species","Classification","TE_Presence","Gene_Count")
	for(x in unique(df7$Classification)){
		df9 <- rbind(df9,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
			Duplicate.2=gsub(".*-","",x),table(df7[df7$Classification == x,c(12,17)])))
	}
	colnames(df9) <- c("Species","Classification","Duplicate.1","Duplicate.2"
						,"Duplicate.1_TE_presence","Duplicate.2_TE_presence","Gene_Count")
	df10 <- df9[df9$Species==i,]
	#df10 will have all the 1-0 genes for gbM-gbM, unM-unM, and teM-teM set to 0-1
	df10[df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==0 & 
		df10$Duplicate.2_TE_presence==1,]$Gene_Count <- df10[df10$Duplicate.1==df10$Duplicate.2 & 
			df10$Duplicate.1_TE_presence==0 & df10$Duplicate.2_TE_presence==1,]$Gene_Count + 
				df10[df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==1 & 
					df10$Duplicate.2_TE_presence==0,]$Gene_Count
	df10[df10$Duplicate.1==df10$Duplicate.2 & df10$Duplicate.1_TE_presence==1 & 
	df10$Duplicate.2_TE_presence==0,]$Gene_Count <- 0
	df10$Gene_Percent <- 0
	for(x in 1:nrow(df10)){
		df10[x,]$Gene_Percent <- df10[x,]$Gene_Count/sum(df10[df10$Species==df10[x,]$Species &
									df10$Classification==df10[x,]$Classification,]$Gene_Count)
	}
	df10$TE_Group <- paste(df10$Duplicate.1_TE_presence,df10$Duplicate.2_TE_presence,sep="-")
	df11 <- rbind(df11,df10)
	df12 <- merge(df5,df6,by.x="Feature",by.y="Duplicate.1")
	df13 <- merge(df5,df6,by.x="Feature",by.y="Duplicate.2")
	colnames(df13) <- colnames(df12)
	df14 <- rbind(df12,df13)
	df15 <- unique(df14[c(1,2,7,8,12,13,14)])
	df16 <- data.frame(Species=i,table(df15$Classification.x,df15$TE2,df15$Classification.y))
	df16$Gene_Percent <- 0
	for(x in 1:nrow(df16)){
		df16[x,]$Gene_Percent <- df16[x,]$Freq/sum(df16[df16$Var3==df16[x,]$Var3 & df16$Var1==df16[x,]$Var1,]$Freq)
	}
	df16 <- df16[df16$Gene_Percent != "NaN",]
	colnames(df16) <- c("Species","Methylation","TE_Presence","Duplicate_Methylation","Gene_Count","Gene_Percent")
	df17 <- rbind(df17,df16)
	#Plot df17
	for(x in c("gbM","unM","teM")){
		#Make the plot
		p <- ggplot(df16[df16$Methylation=="unM",]) + 
			geom_bar(aes(x=Duplicate_Methylation,y=Gene_Percent,fill=TE_Presence),stat="identity",position="dodge")
		#Save the plot
		ggsave(paste(path3,i,"_",x,"_TE_presence.pdf",sep=""),p)
	} #End loop
}





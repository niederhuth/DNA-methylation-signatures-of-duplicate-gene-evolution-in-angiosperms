library(scales)
library(ggplot2)
library(stats)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

species = c("Aduranensis", "Aipaensis", "Alyrata","Athaliana","Atrichopoda",
            "Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
            "Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
            "Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
            "Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
            "Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
            "Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
            "Vvinifera","Zmays")

df17 <- df11 <- data.frame()
for(i in species){
	path1 <- paste(i,"/methylpy/results/",i,sep="")
	path2 <- paste(i,"/dupgen/results-unique/",sep="")
	path3 <- paste("../figures_tables/",i,"/",sep="")
	df1 <- read.table(paste(path1,"_classified_genes.tsv",sep=""),
		header=TRUE,sep="\t")[c("Feature","Classification")]
	df1$Classification <- gsub("Unmethylated","unM",df1$Classification)
	df1 <- df1[df1$Classification %in% c("gbM","teM","unM"),]
	df2 <- read.table(paste(path2,"classified_genes.tsv",sep=""),header=TRUE,sep="\t")
	df2 <- df2[df2$Duplication != "unclassified" & df2$Duplication != "singletons",]
	df3 <- read.table(paste(path1,"_TE_intersections.tsv",sep=""),header=TRUE,sep="\t")
	#df3$Gene <- gsub("araip.Araip", "araip.K30076.gnm1.ann1.Araip", df3[,1])
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
	df6 <- read.csv(paste(path3,i,"_duplicate_pair_met.csv",sep=""),header=TRUE)
	df6 <- df6[df6$Similarity != "Undetermined",]
	df7 <- merge(df6,df5[c(1,4:8)],by.x="Duplicate.1",by.y="Feature",all.x=TRUE)
	df7 <- merge(df7,df5[c(1,4:8)],by.x="Duplicate.2",by.y="Feature",all.x=TRUE)[c(2,3,1,4:17)]
	df9 <- df8 <- data.frame()
	df8 <- rbind(df8,data.frame(Species=i,data.frame(table(df5$Classification,df5$TE2))))
	colnames(df8) <- c("Species","Classification","TE_Presence","Gene_Count")
	for(x in unique(df7$Classification)){
		if(nrow(data.frame(table(df7[df7$Classification == x,c(12,17)]))) == 4){
			df9 <- rbind(df9,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
				Duplicate.2=gsub(".*-","",x),table(df7[df7$Classification == x,c(12,17)])))
		} else {
			tmp <- data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
					Duplicate.2=gsub(".*-","",x),table(df7[df7$Classification == x,c(12,17)]))
			if(nrow(tmp[tmp$TE2.x == 0 & tmp$TE2.y == 0,]) == 0){
				tmp <- rbind(tmp,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
					Duplicate.2=gsub(".*-","",x),TE2.x=as.factor(0),TE2.y=as.factor(0),Freq=0))
			} 
			if(nrow(tmp[tmp$TE2.x == 1 & tmp$TE2.y == 0,]) == 0){
				tmp <- rbind(tmp,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
					Duplicate.2=gsub(".*-","",x),TE2.x=as.factor(1),TE2.y=as.factor(0),Freq=0))
			}
			if(nrow(tmp[tmp$TE2.x == 0 & tmp$TE2.y == 1,]) == 0){
				tmp <- rbind(tmp,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
					Duplicate.2=gsub(".*-","",x),TE2.x=as.factor(0),TE2.y=as.factor(1),Freq=0))
			}
			if(nrow(tmp[tmp$TE2.x == 0 & tmp$TE2.y == 1,]) == 0){
				tmp <- rbind(tmp,data.frame(Species=i,Classification=x,Duplicate.1=gsub("-.*","",x),
					Duplicate.2=gsub(".*-","",x),TE2.x=as.factor(1),TE2.y=as.factor(1),Freq=0))
			}
			df9 <- rbind(df9,tmp)
		}
	}
	p <- ggplot(df8) +
	  geom_bar(aes(x=Classification,y=Gene_Count,fill=TE_Presence),stat="identity", position = "dodge") +
	  theme_bw() +
	  scale_y_continuous("Number of genes", expand=c(0.01,0.01)) +
	  scale_fill_discrete(name = "Transposons", labels = c("Absent", "Present"))
	q <- ggplot(df8) +
	  geom_bar(aes(x=Classification,y=Gene_Count,fill=TE_Presence),stat="identity", position = "fill") +
	  theme_bw() +
	  scale_y_continuous("Percent genes", expand=c(0.01,0.01), labels = percent) +
	  scale_fill_discrete(name = "Transposons", labels = c("Absent", "Present"))
	ggsave(paste(path3,"/",i,"_TE_intersection_dodge.pdf",sep=""),p,device="pdf", width = 4, height = 3)
	ggsave(paste(path3,"/",i,"_TE_intersection_percent.pdf",sep=""),q,device="pdf", width = 4, height = 3)
	write.csv(df8,paste(path3,i,"_TE_intersection.csv",sep=""),quote=FALSE,row.names=FALSE)
	colnames(df9) <- c("Species","Classification","Duplicate.1","Duplicate.2"
						,"Duplicate.1_TE_presence","Duplicate.2_TE_presence","Gene_Count")
	#df10 will have all the 1-0 genes for gbM-gbM, unM-unM, and teM-teM set to 0-1
	df10 <- df9
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
	write.csv(df10,paste(path3,i,"_TE_intersection_classification.csv",sep=""),quote=FALSE,row.names=FALSE)
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
	write.csv(df16,paste(path3,i,"_TE_intersection_switching.csv",sep=""),quote=FALSE,row.names=FALSE)
	df17 <- rbind(df17,df16)
	#Plot df16
	for(x in c("gbM","unM","teM")){
		#Make the plot
		p <- ggplot(df16[df16$Methylation==x,]) +
		  theme_bw()+
			geom_bar(aes(x=Duplicate_Methylation,y=Gene_Percent,fill=TE_Presence),stat="identity",position="fill")+
		  scale_y_continuous("Percent genes", expand=c(0.01,0.01), labels = percent) +
		  scale_fill_discrete(name = "Transposons", labels = c("Absent", "Present"))
		#Save the plot
		ggsave(paste(path3,i,"_",x,"_TE_presence-fill.pdf",sep=""),p, width = 4, height = 3)
	}
}
#save tables
path4 <- "../figures_tables/TE_intersections/"
if(!file.exists(path4)){
	dir.create(path4)
}
write.csv(df11,paste(path4,"All_TE_intersection_classification.csv",sep=""),quote=FALSE,row.names=FALSE)
write.csv(df17,paste(path4,"All_TE_intersection_switching.csv",sep=""),quote=FALSE,row.names=FALSE)

df17$p.value <- df17$OR <- NA
for(i in c(1:nrow(df17))){
	tmp <- fisher.test(matrix(c(
				sum(df17[df17$Species==df17[i,]$Species & df17$Methylation==df17[i,]$Methylation & df17$TE_Presence==1,]$Gene_Count),
				sum(df17[df17$Species==df17[i,]$Species & df17$Methylation != df17[i,]$Methylation & df17$TE_Presence==1,]$Gene_Count),
				sum(df17[df17$Species==df17[i,]$Species & df17$Methylation==df17[i,]$Methylation & df17$TE_Presence==0,]$Gene_Count),
				sum(df17[df17$Species==df17[i,]$Species & df17$Methylation != df17[i,]$Methylation & df17$TE_Presence==0,]$Gene_Count)
			),c(2,2)),alternative ="two.sided")
	df17[i,]$OR <- tmp$estimate
	df17[i,]$p.value <- tmp$p.value
}

df18 <-data.frame()
for(a in species){
	for(b in c("gbM","unM","teM")){
		df18 <- rbind(df18,data.frame(Species=a,Methylation=b,
			OR=unique(df17[df17$Species==a & df17$Methylation==b,]$OR),
			p.value=unique(df17[df17$Species==a & df17$Methylation==b,]$p.value)))
	}
}
df18$p.adjust <- p.adjust(df18$p.value,method="BH")
write.csv(df18,paste(path4,"TE_enrichment.csv",sep=""),quote=FALSE,row.names=FALSE)

#Make tables for heatmap
#If starting from table, then start from next line
df18 <- read.csv(paste(path4,"TE_enrichment.csv",sep=""),header=TRUE,row.names=1)
df20 <- df19 <-data.frame()
for(a in species){
	df19 <- rbind(df19,data.frame(Species=a,gbM=df18[df18$Species==a & df18$Methylation=="gbM",]$OR,
		unM=df18[df18$Species==a & df18$Methylation=="unM",]$OR,
		teM=df18[df18$Species==a & df18$Methylation=="teM",]$OR))
	df20 <- rbind(df20,data.frame(Species=a,gbM=df18[df18$Species==a & df18$Methylation=="gbM",]$p.adjust,
		unM=df18[df18$Species==a & df18$Methylation=="unM",]$p.adjust,
		teM=df18[df18$Species==a & df18$Methylation=="teM",]$p.adjust))
}
rownames(df19) <- df19$Species
rownames(df20) <- df20$Species
df20$gbM <- ifelse(df20$gbM < 0.05,NA,"NS")
df20$unM <- ifelse(df20$unM < 0.05,NA,"NS")
df20$teM <- ifelse(df20$teM < 0.05,NA,"NS")
df20$order <- df19$order <- c(35,34,20,19,1,5,23,22,40,14,17,26,27,18,28,13,2,21,32,33,37,
	16,38,3,31,24,43,39,4,11,29,25,10,36,30,7,8,41,42,9,15,12,6)
df19 <- arrange(df19,order)
df20 <- arrange(df20,order)
#TE Enrichment Heatmap
pdf(paste(path4,"TE_enrichment.pdf",sep=""))
heatmap.2(as.matrix(df19[,c(2:4)]), 
          dendrogram='none', 
          cellnote=as.matrix(df20[,c(2:4)]), notecol="black", notecex = 0.75,
          Rowv=FALSE, Colv= FALSE, key = FALSE,      
          breaks = c(0, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.00, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2),
          density.info = "none", trace ="none", 
          sepcolor = TRUE,
          srtCol=60, cexRow = 1, cexCol = 1, offsetCol = FALSE,
          col = cm.colors(14), margins=c(6,8))
dev.off()




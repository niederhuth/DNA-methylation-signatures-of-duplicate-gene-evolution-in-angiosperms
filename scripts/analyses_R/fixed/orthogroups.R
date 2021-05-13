library(ggplot2)
library(reshape2)
library(scales)
library(ape)
library(pheatmap)
library(Cairo)

#List of all species
species_list=c("Acoerulea","Aduranensis","Aipaensis","Alyrata","Aofficinalis","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cannuum","Cclementina","Clanatus","Cmelo",
	"Cpapaya","Crubella","Csativus","Cviolacea","Egrandis","Eguineensis","Esalsugineum",
	"Fvesca","Fxananassa","Gmax","Graimondii","Gsoja","Ljaponicus","Macuminata","Mdomestica",
	"Mesculenta","Mguttatus","Mtruncatula","Nattenuata","Nnucifera","Osativa","Othomaeum",
	"Paxillaris","Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Rchinensis","Roccidentalis","Sbicolor","Sitalica","Slycopersicum","Smelongena","Spolyrhiza",
	"Stuberosum","Sviridis","Tcacao","Vcorymbosum","Vvinifera","Zjujuba","Zmays")
#Dataframe of species with methylation data and their order based on phylogenetic tree
species <- data.frame(Species=c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
	"Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
	"Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
	"Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
	"Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
	"Vvinifera","Zmays"),Order=c(35,34,20,19,1,5,23,22,40,14,17,26,27,18,28,13,2,21,32,33,37,
	16,38,3,31,24,43,39,4,11,29,25,10,36,30,7,8,41,42,9,15,12,6))
#Dataframe of angiosperm families with more than 2 representative species
families <- data.frame(Family=c(rep("Poaceae",9),rep("Brassicaceae",6),rep("Fabaceae",7),
	rep("Rosaceae",7),rep("Solanaceae",6),rep("Cucurbitaceae",3)),
	Species=c("Bdistachyon","Osativa","Phallii","Pvirgatum","Sbicolor","Sitalica","Sviridis",
		"Zmays","Othomaeum","Alyrata","Athaliana","Brapa","Boleracea","Crubella","Esalsugineum",
		"Aduranensis","Aipaensis","Gmax","Gsoja","Ljaponicus","Mtruncatula","Pvulgaris","Fvesca",
		"Fxananassa","Mdomestica","Ppersica","Pxbretschneideri","Rchinensis","Roccidentalis",
		"Cannuum","Nattenuata","Paxillaris","Slycopersicum","Smelongena","Stuberosum","Clanatus",
		"Cmelo","Csativus"))

#Read in orthogroup gene counts for all species & reformat
tmp <- setNames(data.frame(matrix(ncol=1,nrow=0)),c("Orthogroup"))
for(a in species_list){
	path1 <- paste(a,"/ref/mcscanx/",a,"_orthogroups.tsv",sep="")
	tmp2 <- read.csv(path1,header=FALSE,sep="\t")
	tmp3 <- data.frame(table(tmp2$V2))
	colnames(tmp3) <- c('Orthogroup',a)
	tmp <- merge(tmp,tmp3,by="Orthogroup",all=TRUE)
}
tmp[is.na(tmp)] <- 0
geneCounts <- tmp[,-1]
rownames(geneCounts) <- tmp[,1]
rm(tmp,tmp2,tmp3,a)
geneCounts$Total <- rowSums(geneCounts)
#Create two new temporary dataframes to identify missing orthogroups and single-copy orthogroups
df1 <- PA <- subset(geneCounts,select = -Total )
#If a species is missing an orthogroup label it as a 0, if the orthogroup is present
#label it as a 1, regardless of how many genes are in that orthogroup
for(i in colnames(PA)){
  PA[i] <- ifelse(geneCounts[i]==0,0,1)
}
#Sum up the row to get how many species have or are missing an orthogroup 
PA$Total <- rowSums(PA)
#If an orthogroup is single-copy in a species, label it as a 1, otherwise label it a 0.
for(i in colnames(df1)){
  df1[i] <- ifelse(geneCounts[i]==1,1,0)
}
#Get the percentage of species where an orthogroup is single-copy
#First sum up the rows
df1$Total <- rowSums(df1)
#We then divide by the PA$Total, to account for species missing an orthogroup
df1$Total <- df1$Total/PA$Total
#Plot out the distribution of species per orthogroup. We will use this to decide what the
#the minimal number of species at which we classify an orthogroup as being part of the core
#angiosperm gene set. This is similar to what was done in Li et al. 2016
#First check for output directory path
path2 <- "../figures_tables/orthogroups"
if(!file.exists(path2)){
	dir.create(path2)
}
#Make the plot & save
p <- ggplot(PA) + 
	geom_histogram(aes(x=Total,fill=cut(Total,c(50,58))),bins=58) + 
	theme_bw() + theme(legend.position="none") + 
	scale_y_continuous("Number of Orthogroups",expand=c(0,0)) + 
	scale_x_continuous("Number of Species in Orthogroup",expand=c(0,0))
ggsave(paste(path2,"/Orthogroup_distribution.pdf",sep=""))
#Make a version without species-specific orthogroups
p <- ggplot(PA[PA$Total > 1,]) + 
	geom_histogram(aes(x=Total,fill=cut(Total,c(50,58))),bins=58) + 
	theme_bw() + theme(legend.position="none") + 
	scale_y_continuous("Number of Orthogroups",expand=c(0,0)) + 
	scale_x_continuous("Number of Species in Orthogroup",expand=c(0,0),
		breaks=c(seq(2,58,8)))
ggsave(paste(path2,"/Orthogroup_distribution_no_species_specific.pdf",sep=""))
#Make a dataframe of "core" angiosperm genes
coreGenes <- geneCounts[row.names(geneCounts) %in% row.names(PA[PA$Total >= 51,]),]
#Temporary dataframe of coregenes for identifying singleCopy genes
df2 <- df1[row.names(df1) %in% row.names(PA[PA$Total >= 51,]),]
#Make a dataframe of single-copy core angiosperm genes. We will allow for a certain percentage 
#of these to be present as multicopy in some species, due to recent WGD and other duplication 
#events
singleCopy <- geneCounts[row.names(geneCounts) %in% 
	row.names(df2[df2$Total >= 0.7,]),]
#
pSC <- data.frame()
for(i in colnames(singleCopy[1:58])){
	pSC <- rbind(pSC,data.frame(Species=i,
		OpSC=data.frame(table(singleCopy[i]))[2,2]/nrow(singleCopy),
		OgbM=NA,OteM=NA,OunM=NA,OUnclassified=NA,OMissing=NA,
		OpgbM=NA,OpteM=NA,OpunM=NA,OpUnclassified=NA,OpMissing=NA,
		gbM=NA,teM=NA,unM=NA,Unclassified=NA,Missing=NA,Total=NA,
		pgbM=NA,pteM=NA,punM=NA,pUnclassified=NA,pMissing=NA))
}
#Identify multiCopy core angiosperm genes
multiCopy <- geneCounts[row.names(geneCounts) %in% row.names(coreGenes) & 
	!(row.names(geneCounts) %in% row.names(singleCopy)),]
#Identify species specific orthogroups
speciesSpecific <- geneCounts[row.names(geneCounts) %in% row.names(PA[PA$Total==1,]),]
#Identify family specific orthogroups
lineageSpecific <- data.frame()
for(i in unique(families$Family)){
	lineageSpecific <- rbind(lineageSpecific,geneCounts[row.names(geneCounts) %in%
		row.names(PA[rowSums(PA[,!(colnames(PA) %in% c(as.character(families[
		families$Family==i,]$Species),"Total"))]) < 1 & PA$Total >= 2,]),])
}
#Identify orthogroups that don't fit any of the above defitions
otherOG <- geneCounts[!(row.names(geneCounts) %in% c(row.names(coreGenes),
	row.names(speciesSpecific),row.names(lineageSpecific))),]
#Create dataframe of each orthogroup's category
ogCat <- data.frame(Orthogroup=row.names(geneCounts),ogCat=NA)
ogCat$ogCat <- ifelse(ogCat$Orthogroup %in% row.names(speciesSpecific),
	"Species/Lineage Specific",ogCat$ogCat)
ogCat$ogCat <- ifelse(ogCat$Orthogroup %in% row.names(lineageSpecific),
	"Family Specific",ogCat$ogCat)
ogCat$ogCat <- ifelse(ogCat$Orthogroup %in% row.names(otherOG),
	"Cross-Family",ogCat$ogCat)
ogCat$ogCat <- ifelse(ogCat$Orthogroup %in% row.names(multiCopy),
	"Core: Other",ogCat$ogCat)
ogCat$ogCat <- ifelse(ogCat$Orthogroup %in% row.names(singleCopy),
	"Core: Single Copy",ogCat$ogCat)

top <- pOG <- data.frame()
#Iterate over each species
for(a in species$Species){
	#Read in the list of genes classified by methylation status
	df1 <- read.table(paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep=""),
		sep="\t",header=TRUE)
	#Change Unmethylated to unM
	df1$Classification <- gsub("Unmethylated","unM",df1$Classification)
	#Set classification to character, so that R doesn't screw it up
	df1$Classification <- as.character(df1$Classification)
	#Change "NA" to "Missing"
	df1$Classification <- ifelse(is.na(df1$Classification),"Missing",df1$Classification)
	#Read in the orthogroups for that species
	df2 <- read.table(paste(a,"/ref/mcscanx/",a,"_orthogroups.tsv",sep=""),
		header=FALSE,sep="\t")[c(1,2)]
	#Redo the original orthogroup list for that species for other analyses
	new <- merge(df2,ogCat,by.x="V2",by.y="Orthogroup")[c("V1","V2","ogCat")]
	write.table(new,paste(a,"/ref/mcscanx/",a,"_orthogroups.tsv",sep=""),
		col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
	#Merge DNA methylation classification & orthogroups
	df3 <- merge(df2,df1[c(1,30)],by.x="V1",by.y="Feature")
	#Table the values so that we can get the count of each classification for each orthogroup
	df4 <- data.frame(table(df3[c(2,3)]))
	#Build a new table with each column being the DNA methylation classification
	#and each row an orthogroup
	df5 <- df4[df4$Classification=="gbM",c(1,3)]
	df5 <- merge(df5,df4[df4$Classification=="teM",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="unM",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Unclassified",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Missing",c(1,3)],by="V2")
	#Rename columns
	colnames(df5) <- c("Orthogroup","gbM","teM","unM","Unclassified","Missing")
	#Sum up each row
	df5$Total <- rowSums(df5[c(2:6)])
	#Reformat for easier plotting
	df6 <- melt(df5)
	#Populate the pSC dataframe for analysis
	tmp <- data.frame(t(colSums(df5[df5$Orthogroup %in% row.names(singleCopy),c(2:7)])))
	for(b in c('gbM','teM','unM','Unclassified','Missing','Total')){
		pSC[pSC$Species==a,b] <- tmp[b]
		if(b != 'Total'){
			pSC[pSC$Species==a,paste("p",b,sep="")] <- tmp[b]/tmp['Total']
			pSC[pSC$Species==a,paste("O",b,sep="")] <- 
			nrow(df5[df5$Orthogroup %in% row.names(singleCopy) & df5[b] > 0,])
			pSC[pSC$Species==a,paste("Op",b,sep="")] <- 
			nrow(df5[df5$Orthogroup %in% row.names(singleCopy) & df5[b] > 0,])/
			nrow(singleCopy)
			#Make a dataframe of the top 5 orthogroups with the most of each 
			#methylation classification
			top <- rbind(top,data.frame(Species=a,Classification=b,
				Orthogroup=df5[df5[[b]] %in% sort(df5[[b]],decreasing=T)[1:5],"Orthogroup"],
				Count=df5[df5[[b]] %in% sort(df5[[b]],decreasing=T)[1:5],b]))
		}
	}
	rm(tmp)
	#Create a dataframe counting genes by methylation class and different categories of
	#orthogroups
	df7 <- merge(ogCat,df6,by="Orthogroup")
	df7$Species <- a
	df7$Order <- species[species$Species==a,]$Order
	df7$pOG <- df7$pmC <- NA
	for(x in unique(df7$ogCat)){
		for(y in unique(df7$variable)){
			df7[df7$ogCat==x & df7$variable==y,]$pmC <- 
				df7[df7$ogCat==x & df7$variable==y,]$value/sum(df7[df7$variable==y,]$value)
			df7[df7$ogCat==x & df7$variable==y,]$pOG <- 
				df7[df7$ogCat==x & df7$variable==y,]$value/
				sum(df7[df7$ogCat==x & df7$variable != "Total",]$value)
		}
	} 
	pOG <- rbind(pOG,df7)
}

#Top orthogroups
top <- merge(top,ogCat,by="Orthogroup")
write.csv(top,paste(path2,"/","top_orthogroups.csv",sep=""),quote=FALSE,row.names=FALSE)

#This is an ugly workaround to reduce the size of the graphs by summarizing data
pOG2 <- data.frame(Species=c(),Order=c(),ogCat=c(),mC=c(),Count=c(),pOG=c(),pmC=c())
for(a in species$Species){
	for(b in unique(pOG$ogCat)){
		for(c in unique(pOG$variable)){
			pOG2 <- rbind(pOG2,data.frame(Species=a,Order=species[species$Species == a,]$Order,
				ogCat=b,mC=c,Count=sum(pOG[pOG$Species == a & pOG$ogCat == b & pOG$variable == c,]$value),
				pOG=sum(pOG[pOG$Species == a & pOG$ogCat == b & pOG$variable == c,]$value)/
					sum(pOG[pOG$Species == a & pOG$ogCat == b & pOG$variable == "Total",]$value),
				pmC=sum(pOG[pOG$Species == a & pOG$ogCat == b & pOG$variable == c,]$value)/
					sum(pOG[pOG$Species == a & pOG$variable == c,]$value)))
		}
	}
}
pOG2$Order2 <- c(rep(5,6),rep(3,6),rep(1,6),rep(2,6),rep(4,6))
pOG2$Order3 <- c(6,3,2,5,1)

#Fisher's Exact test for enrichment/depletion
pOG2$p.value <- pOG2$OR <- NA
for(i in 1:nrow(pOG2)){
	if(!(pOG2[i,]$mC %in% c("Total","Missing","Unclassified"))){
		tmp <- fisher.test(matrix(c(
			pOG2[pOG2$Species==pOG2[i,]$Species & pOG2$ogCat==pOG2[i,]$ogCat & pOG2$mC==pOG2[i,]$mC,]$Count,
			pOG2[pOG2$Species==pOG2[i,]$Species & pOG2$ogCat==pOG2[i,]$ogCat & pOG2$mC == "Total",]$Count -
			pOG2[pOG2$Species==pOG2[i,]$Species & pOG2$ogCat==pOG2[i,]$ogCat & pOG2$mC==pOG2[i,]$mC,]$Count,
			sum(pOG2[pOG2$Species==pOG2[i,]$Species & pOG2$mC==pOG2[i,]$mC & pOG2$ogCat!=pOG2[i,]$ogCat,]$Count),
			sum(pOG2[pOG2$Species==pOG2[i,]$Species & pOG2$ogCat!=pOG2[i,]$ogCat & 
				pOG2$mC!=pOG2[i,]$mC & pOG2$mC != "Total",]$Count)
			),c(2,2)),alternative ="two.sided")
		pOG2[i,]$OR <- tmp$estimate
		pOG2[i,]$p.value <- tmp$p.value
		} else {
			pOG2[i,]$OR <- NA
			pOG2[i,]$p.value <- NA
		}
}
rm(tmp)
pOG2$p.adjust <- p.adjust(pOG2$p.value,method="BH")
pOG3 <- na.omit(pOG2[,c(1,3,4,5,10,11)])
write.csv(pOG3,paste(path2,"/enrichment_test.csv",sep=""),row.names=FALSE,quote=FALSE)

#Plot percentage of each methylation class that is in each type of orthogroup
for(a in c("gbM","teM","unM","Unclassified","Missing")){
	p <- ggplot(pOG2[pOG2$mC==a,]) + 
		geom_bar(aes(x=reorder(Species,Order),y=pmC,fill=reorder(ogCat,Order2)),stat="identity") + 
		theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
			axis.title.x=element_blank()) + 
		scale_y_continuous("Percentage of Genes",expand=c(0,0),labels=percent) +
		scale_fill_discrete(name = "Orthogroup Class")
	ggsave(paste(path2,"/",a,"_orthogroup_distribution.pdf",sep=""),p,width=10,height=4)
}
#Plot percentage of each type of orthogroup that is in each methylation class
for(a in c("Species/Lineage Specific","Family Specific","Cross-Family","Core: Other",
	"Core: Single Copy")){
	p <- ggplot(pOG2[pOG2$ogCat==a & pOG2$mC != "Total",]) + 
		geom_bar(aes(x=reorder(Species,Order),y=pOG,fill=reorder(mC,Order3)),stat="identity") + 
		theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
			axis.title.x=element_blank()) + 
		scale_y_continuous("Percentage of Genes",expand=c(0,0),labels=percent) +
		scale_fill_discrete(name = "Methylation Class")
	ggsave(paste(path2,"/",gsub(":","",gsub("/","-",gsub(" ","_",a))),"_methylation_distribution.pdf",
		sep=""),p,width=10,height=4)
}

#Examine core gene clustering based on methylation state
cm <- data.frame(Orthogroup=unique(pOG[pOG$ogCat=="Core: Single Copy" | pOG$ogCat=="Core: Other",]$Orthogroup))
for(a in c("gbM","teM","unM","Unclassified","Missing")){
	tmp <- data.frame(table(as.vector(pOG[pOG$ogCat=="Core: Single Copy" & pOG$variable==a 
		& pOG$value > 0 | pOG$ogCat=="Core: Other" & pOG$variable==a & pOG$value > 0,]$Orthogroup)))
	colnames(tmp) <- c("Orthogroup",a)
	cm <- merge(cm,tmp,by="Orthogroup",all=TRUE)
}
rm(tmp)
cm[is.na(cm)] <- 0
cm$Total <- rowSums(cm[c(2:6)])
cm[8:12] <- cm[2:6]
colnames(cm) <- c("Orthogroup","gbM","teM","unM","Unclassified","Missing","Total",
	"percent_gbM","percent_teM","percent_unM","percent_Unclassified","percent_Missing")
cm[c(8:12)] <- cm[c(8:12)]/cm$Total
cm2 <- melt(cm[c(8:12)])
cmHTM <- pheatmap(cm[c(8:12)])
cairo_pdf(paste(path2,'/core_orthogroup_cluster_heatmap.pdf',sep=''),family="Arial")
	pheatmap(cm[c(8:12)],cutree_rows=3,show_rownames=FALSE,
		legend_breaks=c(0.2,0.4,0.6,0.8),legend_labels=c("20%","40%","60%","80%"),
		labels_col=c("gbM","teM","unM","Unclassified","Missing"))
dev.off()
cm <- cbind(cm,cluster=cutree(cmHTM$tree_row,3))
write.table(cm,paste(path2,"/core_orthogroups.tsv",sep=""),quote=FALSE,row.names=FALSE,sep="\t")

#Examine single copy gene clustering based on methylation state
scm <- data.frame(Orthogroup=unique(pOG[pOG$ogCat=="Core: Single Copy",]$Orthogroup))
for(a in c("gbM","teM","unM","Unclassified","Missing")){
	tmp <- data.frame(table(as.vector(pOG[pOG$ogCat=="Core: Single Copy" & 
		pOG$variable==a & pOG$value > 0,]$Orthogroup)))
	colnames(tmp) <- c("Orthogroup",a)
	scm <- merge(scm,tmp,by="Orthogroup",all=TRUE)
}
rm(tmp)
scm[is.na(scm)] <- 0
scm$Total <- rowSums(scm[c(2:6)])
scm[8:12] <- scm[2:6]
colnames(scm) <- c("Orthogroup","gbM","teM","unM","Unclassified","Missing","Total",
	"percent_gbM","percent_teM","percent_unM","percent_Unclassified","percent_Missing")
scm[c(8:12)] <- scm[c(8:12)]/cm$Total
scmHTM <- pheatmap(scm[c(8:12)])
cairo_pdf(paste(path2,'/single_copy_cluster_heatmap.pdf',sep=''),family="Arial")
	pheatmap(scm[c(8:12)],cutree_rows=3,show_rownames=FALSE,
		legend_breaks=c(0.2,0.4,0.6,0.8),legend_labels=c("20%","40%","60%","80%"),
		labels_col=c("gbM","teM","unM","Unclassified","Missing"))
dev.off()
scm <- cbind(scm,cluster=cutree(scmHTM$tree_row,3))
write.table(scm,paste(path2,"/single_copy_orthogroups.tsv",sep=""),quote=FALSE,row.names=FALSE,sep="\t")

#Run GO term enrichment on the different clusters of core and single copy genes, using Arabidopsis
library(topGO)
library(GO.db)
source("../scripts/analyses_R/functions.R")
path3 <- "../figures_tables/orthogroups/GO_terms/"
if(!file.exists(path3)){
	dir.create(path3)
}
#Read in GO terms for Arabidopsis in topGO format
goTerms <- readMappings(file="Athaliana/ref/annotations/Athaliana-topGO.txt")
#Read in orthogroup-gene mappings for Arabidopsis
Og2At <- read.table("Athaliana/ref/mcscanx/Athaliana_orthogroups.tsv",header=FALSE,sep="\t")
#Core unM
cm2AtUM <- factor(as.integer(Og2At$V1 %in% Og2At[Og2At$V2 %in% cm[cm$cluster==3,]$Orthogroup,]$V1))
names(cm2AtUM) <- Og2At$V1
topGO(cm2AtUM,goTerms,nodeSize=5,fdr=0.05,filename="Core_UM",path=path3,returnData=FALSE)
#Core gbM
cm2AtgBM <- factor(as.integer(Og2At$V1 %in% Og2At[Og2At$V2 %in% cm[cm$cluster==2,]$Orthogroup,]$V1))
names(cm2AtgBM) <- Og2At$V1
topGO(cm2AtgBM,goTerms,nodeSize=5,fdr=0.05,filename="Core_gbM",path=path3,returnData=FALSE)
#Core Mixed methylation
cm2AtMixed <- factor(as.integer(Og2At$V1 %in% Og2At[Og2At$V2 %in% cm[cm$cluster==1,]$Orthogroup,]$V1))
names(cm2AtMixed) <- Og2At$V1
topGO(cm2AtMixed,goTerms,nodeSize=5,fdr=0.05,filename="Core_Mixed",path=path3,returnData=FALSE)
#Core Single Copy unM
scm2AtUM <- factor(as.integer(Og2At$V1 %in% Og2At[Og2At$V2 %in% scm[scm$cluster==3,]$Orthogroup,]$V1))
names(scm2AtUM) <- Og2At$V1
topGO(scm2AtUM,goTerms,nodeSize=5,fdr=0.05,filename="single_copy_UM",path=path3,returnData=FALSE)
#Core Single Copy gBM
scm2AtgBM <- factor(as.integer(Og2At$V1 %in% Og2At[Og2At$V2 %in% scm[scm$cluster==2,]$Orthogroup,]$V1))
names(scm2AtgBM) <- Og2At$V1
topGO(scm2AtgBM,goTerms,nodeSize=5,fdr=0.05,filename="single_copy_gbM",path=path3,returnData=FALSE)
#Core Single Copy Mixed methylation
scm2AtMixed <- factor(as.integer(Og2At$V1 %in% Og2At[Og2At$V2 %in% scm[scm$cluster==1,]$Orthogroup,]$V1))
names(scm2AtMixed) <- Og2At$V1
topGO(scm2AtMixed,goTerms,nodeSize=5,fdr=0.05,filename="single_copy_Mixed",path=path3,returnData=FALSE)


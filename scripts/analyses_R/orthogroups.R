library(ggplot2)
library(reshape2)
library(scales)
library(ape)

#List of species
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

#species=data.frame(Species=c("Athaliana"))

#Set path to orthofinder results
path1=paste("orthofinder/orthofinder/",dir("orthofinder/orthofinder/"),sep="")
#Read in the rooted species tree
tre <- read.tree(paste(path1,"/Species_Tree/SpeciesTree_rooted_node_labels.txt",sep=""))
#Read in orthogroup gene counts for all species & reformat
geneCounts <- read.table(paste(path1,"/Orthogroups/Orthogroups.GeneCount.tsv",sep=""),
	header=T,sep="\t",row.names=1)
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
df1$Total <- rowSums(df1)/ncol(df1)
#Plot out the distribution of species per orthogroup. We will use this to decide what the
#the minimal number of species at which we classify an orthogroup as being part of the core
#angiosperm gene set. This is similar to what was done in Li et al. 2016
p <- ggplot(PA) + 
	geom_histogram(aes(x=Total,fill=cut(Total,c(50,58))),bins=58) + 
	theme_bw() + theme(legend.position="none") + 
	scale_y_continuous("Number of Orthogroups",expand=c(0,0)) + 
	scale_x_continuous("Number of Species in Orthogroup",expand=c(0,0))
#Make a dataframe of "core" angiosperm genes
coreGenes <- geneCounts[row.names(geneCounts) %in% row.names(PA[PA$Total >= 51,]),]
#Temporary dataframe of coregenes for identifying singleCopy genes
df2 <- df1[row.names(df1) %in% row.names(PA[PA$Total >= 51,]),]
#Make a dataframe of single-copy core angiosperm genes. We will allow for a certain percentage 
#of these to be present as multicopy in some species, due to recent WGD and other duplication 
#events
singleCopy <- geneCounts[row.names(geneCounts) %in% 
	row.names(df2[df2$Total >= 0.75,]),]
#
pSC <- data.frame()
for(i in colnames(singleCopy[1:58])){
	pSC <- rbind(pSC,data.frame(Species=i,
		OpSC=data.frame(table(singleCopy[i]))[2,2]/nrow(singleCopy),
		OgbM=NA,OTE.like=NA,OUnmethylated=NA,OUnclassified=NA,OMissing=NA,
		OpgbM=NA,OpTE.like=NA,OpUnmethylated=NA,OpUnclassified=NA,OpMissing=NA,
		gbM=NA,TE.like=NA,Unmethylated=NA,Unclassified=NA,Missing=NA,Total=NA,
		pgbM=NA,pTE.like=NA,pUnmethylated=NA,pUnclassified=NA,pMissing=NA))
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
	df3 <- merge(df2,df1[c(1,23)],by.x="V1",by.y="Feature")
	#Table the values so that we can get the count of each classification for each orthogroup
	df4 <- data.frame(table(df3[c(2,3)]))
	#Build a new table with each column being the DNA methylation classification
	#and each row an orthogroup
	df5 <- df4[df4$Classification=="gbM",c(1,3)]
	df5 <- merge(df5,df4[df4$Classification=="TE-like",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Unmethylated",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Unclassified",c(1,3)],by="V2")
	df5 <- merge(df5,df4[df4$Classification=="Missing",c(1,3)],by="V2")
	#Rename columns
	colnames(df5) <- c("Orthogroup","gbM","TE.like","Unmethylated","Unclassified","Missing")
	#Sum up each row
	df5$Total <- rowSums(df5[c(2:6)])
	#Reformat for easier plotting
	df6 <- melt(df5)
	#Populate the pSC dataframe for analysis
	tmp <- data.frame(t(colSums(df5[df5$Orthogroup %in% row.names(singleCopy),c(2:7)])))
	for(b in c('gbM','TE.like','Unmethylated','Unclassified','Missing','Total')){
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
write.csv(top,"top_orthogroups.csv",quote=FALSE,row.names=FALSE)

#Plot
for(a in c("gbM","TE.like","Unmethylated","Unclassified","Missing")){
	p <- ggplot(pOG[pOG$variable==a,]) + 
		geom_bar(aes(x=reorder(Species,Order),y=pmC,fill=ogCat),stat="identity") + 
		theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
		scale_y_continuous("Percentage of Genes",expand=c(0,0))
	ggsave(paste(a,"_orthogroup_distribution.pdf",sep=""),p,width=10,height=4)
}

for(a in c("Species/Lineage Specific","Family Specific","Cross-Family","Core: Other",
	"Core: Single Copy")){
	p <- ggplot(pOG[pOG$ogCat==a & pOG$variabl != "Total",]) + 
		geom_bar(aes(x=reorder(Species,Order),y=pOG,fill=variable),stat="identity") + 
		theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
		scale_y_continuous("Percentage of Genes",expand=c(0,0))
	ggsave(paste(gsub(":","",gsub("/","-",gsub(" ","_",a))),"_methylation_distribution.pdf",
		sep=""),p,width=10,height=4)
}


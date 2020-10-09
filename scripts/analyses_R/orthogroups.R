library(ggplot2)
library(reshape2)
library(scales)
library(ape)

#List of species
species = c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
	"Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
	"Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
	"Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
	"Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
	"Vvinifera","Zmays")

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

pSC <- data.frame()
for(i in colnames(singleCopy[1:58])){
	pSC <- rbind(pSC,data.frame(Species=i,
		OpSC=data.frame(table(singleCopy[i]))[2,2]/nrow(singleCopy),
		OgbM=NA,OTE.like=NA,OUnmethylated=NA,OUnclassified=NA,OMissing=NA,
		OpgbM=NA,OpTE.like=NA,OpUnmethylated=NA,OpUnclassified=NA,OpMissing=NA,
		gbM=NA,TE.like=NA,Unmethylated=NA,Unclassified=NA,Missing=NA,Total=NA,
		pgbM=NA,pTE.like=NA,pUnmethylated=NA,pUnclassified=NA,pMissing=NA))
}

#Identify species specific orthogroups
speciesSpecific <- geneCounts[row.names(geneCounts) %in% row.names(PA[PA$Total==1,]),]
#Identify high copy number species specific orthogroups. These are suspect of being TEs
highCopySS <- speciesSpecific[speciesSpecific$Total > 20,]

top <- df7 <- data.frame()
#Iterate over each species
for(a in species){
	#Read in the list of genes classified by methylation status
	df1 <- read.table(paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep=""),
		sep="\t",header=TRUE)
	#Set classification to character, so that R doesn't screw it up
	df1$Classification <- as.character(df1$Classification)
	#Change "NA" to "Missing"
	df1$Classification <- ifelse(is.na(df1$Classification),"Missing",df1$Classification)
	#Read in the orthogroups for that species
	df2 <- read.table(paste(a,"/ref/mcscanx/",a,"_orthogroups.tsv",sep=""),header=FALSE,sep="\t")
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
			top <- rbind(top,data.frame(Species=a,Classification=b,
				Orthogroup=df5[df5[[b]] %in% sort(df5[[b]],decreasing=T)[1:5],"Orthogroup"],
				Count=df5[df5[[b]] %in% sort(df5[[b]],decreasing=T)[1:5],b]))
		}
	}
	rm(tmp)
	for(b in c('gbM','TE.like','Unmethylated')){
		if(nrow(df5[df5[b] > 0 & df5$Orthogroup %in% row.names(singleCopy),]) != 0){
			df7 <- rbind(df7,data.frame(Species=c(a),Classification=c(b),
			Count=df5[df5[b] > 0 & df5$Orthogroup %in% row.names(singleCopy),]$Total))
		}
	}
}

#Top orthogroups
write.csv(top,"top_orthogroups.csv",quote=FALSE,row.names=FALSE)





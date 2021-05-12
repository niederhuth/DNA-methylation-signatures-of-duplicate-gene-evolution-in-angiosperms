library(dplyr)
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

#species=c("Athaliana")

#Create empty dataframe for K-S test results
KaKs <- Ks <- data.frame()

#Loop over each species
for(a in species){
	#output path
	path1 <- paste("../figures_tables/",a,sep="")
	if(!file.exists(path1)){
		dir.create(path1)
	}
	#Genes classified by mode of duplication
	path2 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
	#Genes classified by methylation status
	path3 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
	#Read in the files
	df1 <- read.table(path2,header=TRUE,sep="\t")
	df2 <- read.table(path3,header=TRUE,sep="\t")
	#Merge the two into a new dataframe
	df3 <- merge(df1,df2[,c(1,30)],by.x="Feature",by.y="Feature")
	#Read in duplicate similarity
	path4 <- paste("../figures_tables/", a,"/",a, "_Duplicate_pair_met2.csv",sep="")
	df4 <- read.csv(path4,header=TRUE)
	#This ugly addition was added later to combine analyses of all SDGs
	df9 <- data.frame()
	#Loop through each category of gene duplication
	for(b in c("wgd","proximal","dispersed","tandem","transposed")){
		#Read in appropriate kaks results
		path5a <- paste(a,"/dupgen/results/kaks_results","/",a,".",b,".kaks",sep="")
		path5u <- paste(a,"/dupgen/results-unique/kaks_results","/",a,".",b,".kaks",sep="")
		df5a <- read.table(path5a,header=T,sep="\t")
		df5u <- read.table(path5u,header=T,sep="\t")
		df6 <- merge(df4,df5u,by=c("Duplicate.1","Duplicate.2"))
		df7 <- setdiff(df4[df4$Duplication==b,],df6[c(1:7)])
		df8 <- rbind(df6,merge(df7,df5a,by=c("Duplicate.1","Duplicate.2")))




		
		#Plot the density of genes for each methylation class based Ks
		p <- ggplot(df8) +
			geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
			geom_density(aes(x=Ks,color='Total'),size=1) +
			theme_bw() +
			theme(axis.text=element_text(color="black"),
				axis.ticks=element_line(color="black"))+
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
			scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
				"teM"="#E69F00","Unmethylated"="#CC79A7","Unclassified"=NA),
				breaks=c("Total","gbM","teM","Unmethylated"))
		ggsave(paste(path1,"/",a,"_",b,"_Ks.pdf",sep=""),p,device="pdf")
		#Plot the density of genes for each methylation class based Ka/Ks ratio
		p <- ggplot(df8) +
			geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
			geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
			geom_vline(xintercept=1, linetype="longdash", color="grey55") +
			theme_bw() +
			theme(axis.text=element_text(color="black"),
				axis.ticks=element_line(color="black"))+
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
			scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
				"teM"="#E69F00","Unmethylated"="#CC79A7","Unclassified"=NA),
				breaks=c("Total","gbM","teM","Unmethylated"))
		ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
		#Perform K-S test on each distribution
		for(c in c("teM","gbM","Unmethylated")){
			if(nrow(df8[df8$Classification == c,]) != 0){
				tmp <- ks.test(df8[df8$Classification == c,]$Ks,
					sample(df8$Ks,size=nrow(df8[df8$Classification == c,]),replace=FALSE))
				Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
					D.statistic=tmp$stat,p.value=tmp$p.value))
				tmp <- ks.test(df8[df8$Classification == c,]$Ka.Ks,
					sample(df8$Ka.Ks,size=nrow(df8[df8$Classification == c,]),replace=FALSE))
				KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
					D.statistic=tmp$stat,p.value=tmp$p.value))
			} else {
				Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
					D.statistic=NA,p.value=NA))
				KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
					D.statistic=NA,p.value=NA))
			}
		}
		df9 <- rbind(df9,df8)
	}
	#Plot the density of genes for each methylation class based Ks for SDGs
	p <- ggplot(df9[df9$Duplication != "wgd",]) +
		geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
		geom_density(aes(x=Ks,color='Total'),size=1) +
		theme_bw() +
		theme(axis.text=element_text(color="black"),
			axis.ticks=element_line(color="black"))+
		scale_y_continuous(expand=c(0,0)) +
		scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
		scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
			"teM"="#E69F00","Unmethylated"="#CC79A7","Unclassified"=NA),
			breaks=c("Total","gbM","teM","Unmethylated"))
	ggsave(paste(path1,"/",a,"_SDGs_Ks.pdf",sep=""),p,device="pdf")
	#Plot the density of genes for each methylation class based Ka/Ks ratio for SDGs
	p <- ggplot(df9[df9$Duplication != "wgd",]) +
		geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
		geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
		geom_vline(xintercept=1, linetype="longdash", color="grey55") +
		theme_bw() +
		theme(axis.text=element_text(color="black"),
			axis.ticks=element_line(color="black"))+
		scale_y_continuous(expand=c(0,0)) +
		scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
		scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
			"teM"="#E69F00","Unmethylated"="#CC79A7","Unclassified"=NA),
			breaks=c("Total","gbM","teM","Unmethylated"))
	ggsave(paste(path1,"/",a,"_SDGs_KaKs.pdf",sep=""),p,device="pdf")
	#Output csv of Ka/Ks data for each duplicate pair
	write.csv(df9,paste(path1,"/",a,"_KaKs_values.csv",sep=""),quote=FALSE)
}

#Adjust p.values for K-S tests
Ks$p.adjust <- p.adjust(Ks$p.value,method="BH")
KaKs$p.adjust <- p.adjust(KaKs$p.value,method="BH")
#Output K-S results
write.csv(Ks,"../figures_tables/Ks-distribution-test.csv",quote=FALSE)
write.csv(KaKs,"../figures_tables/KaKs-distribution-test.csv",quote=FALSE)











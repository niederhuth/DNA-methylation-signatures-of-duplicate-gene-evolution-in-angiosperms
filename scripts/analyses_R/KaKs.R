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

#Create empty dataframe for K-S test results
KaKs <- Ks <- data.frame()

#Loop over each species
for(a in species){
	#output path
	path1 <- paste("../../figures_tables/",a,sep="")
	#Genes classified by mode of duplication
	path2 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
	#Genes classified by methylation status
	path3 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
	#Read in the files
	df1 <- read.table(path2,header=TRUE,sep="\t")
	df2 <- read.table(path3,header=TRUE,sep="\t")
	#Merge the two into a new dataframe
	df3 <- merge(df1,df2[,c(1,23)],by.x="Feature",by.y="Feature")
	#Loop through each category of gene duplication
	for(b in c("wgd","proximal","dispersed","tandem","transposed")){
		#Read in appropriate kaks results
		path4 <- paste(a,"/dupgen/results/kaks_results","/",a,".",b,".kaks",sep="")
		df4 <- read.table(path4,header=T,sep="\t")
		#Rename colnames so "Duplicate.1" is "Feature"
		colnames(df4) <- c("Feature","Duplicate.2","Ka","Ks","Ka.Ks","P.Value")
		#Because df3 has 1 gene per line, while df4 is for each gene pair,need to 
		#create a duplicate copy of df4 (df5), with columns 1 & 2 inverted and 
		#combine this with df4 so that can merge with df3 without loss of genes
		df5 <- df4[c(2,1,3:6)]
		#Make sure column names match for rbind to properly work
		colnames(df5) <- colnames(df4)
		df6 <- rbind(df4,df5)
		#sometimes, a gene may be represented more than once, because it is ancestor
		#to more than one gene. This will give it multiple Ka & Ks values and result
		#in it being counted more than once. Here we will randomly select one of the
		#values, so as to prevent it from being double counted. Alternatively you can
		#select to use the lowest Ks or highest Ks, just unhash that line and add a 
		#hash to the other.
		df7 <- data.frame()
		for(i in df3[df3$Duplication == b,]$Feature){
			if(nrow(df6[df6$Feature==i,]) != 0){
				df7 <- rbind(df7,
					df6[row.names(df6) == sample(row.names(df6[df6$Feature==i,]),1),])
				#df7 <- rbind(df7,
				#	df6[df6$Feature==i & df6$Ks == min(df6[df6$Feature==i,]$Ks),])
				#df7 <- rbind(df7,
				#	df6[df6$Feature==i & df6$Ks == max(df6[df6$Feature==i,]$Ks),])
			} else {
				df7 <- rbind(df7,data.frame(Feature=i,Duplicate.2=NA,Ka=NA,Ks=NA,
					Ka.Ks=NA,P.Value=NA))
			}
		}
		#Merge with classified genes 
		df8 <- merge(df3,df7,by="Feature")
		#Because some genes in a duplicate pair may have arisen by other mechanisms,
		#need to remove genes in df8 that are not themselves classified as the type of 
		#duplication being analyzed
		df8 <- df8[df8$Duplication==b,]
		#Plot the density of genes for each methylation class based Ks
		p <- ggplot(df8) +
			geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
			geom_density(aes(x=Ks,color='Total'),size=1) +
			theme_bw() +
			theme(axis.text=element_text(color="black"),
				axis.ticks=element_line(color="black"))+
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous("Ks") +
			scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7","Unclassified"="#999999"),
				breaks=c("Total","gbM","TE-like","Unmethylated","Unclassified"))
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
			scale_x_continuous("Ka/Ks") +
			scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7","Unclassified"="#999999"),
				breaks=c("Total","gbM","TE-like","Unmethylated","Unclassified"))
		ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
		#Perform K-S test on each distribution
		for(c in c("TE-like","gbM","Unmethylated")){
			tmp <- ks.test(df8[df8$Classification == c,]$Ks,
				sample(df8$Ks,size=nrow(df8[df8$Classification == c,]),replace=FALSE))
			Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
				D.statistic=tmp$stat,p.value=tmp$p.value))
			tmp <- ks.test(df8[df8$Classification == c,]$Ka.Ks,
				sample(df8$Ka.Ks,size=nrow(df8[df8$Classification == c,]),
					replace=FALSE))
			KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
				D.statistic=tmp$stat,p.value=tmp$p.value))
		}
	}
}

#Adjust p.values for K-S tests
Ks$p.adjust <- p.adjust(Ks$p.value,method="BH")
KaKs$p.adjust <- p.adjust(KaKs$p.value,method="BH")
#Output K-S results
write.csv(Ks,"../../figures_tables/Ks-distribution-test.csv",quote=FALSE)
write.csv(KaKs,"../../figures_tables/KaKs-distribution-test.csv",quote=FALSE)











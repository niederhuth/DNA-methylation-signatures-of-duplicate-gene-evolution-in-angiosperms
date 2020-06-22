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
		path4 <- paste(a,"/dupgen/results-unique/kaks_results","/",a,".",b,".kaks",
			sep="")
		df4 <- read.table(path4,header=T,sep="\t")
		#Because df3 has 1 gene per line, while df4 is for each gene pair,need to 
		#create a duplicate copy of df4 (df5), with columns 1 & 2 inverted and 
		#combine this with df4 so that can merge with df3 without loss of genes
		df5 <- df4[c(2,1,3:6)]
		#Make sure column names match for rbind to properly work
		colnames(df5) <- colnames(df4)
		df6 <- rbind(df4,df5)
		#Merge with classified genes 
		df7 <- merge(df3,df6,by.x="Feature",by.y="Duplicate.1")
		#Because some genes in a duplicate pair may have arisen by other mechanisms,
		#need to remove genes in df7 that are not themselves classified as the type of 
		#duplication being analyzed
		df7 <- df7[df7$Duplication==b,]
		#Plot the density of genes for each methylation class based Ks
		p <- ggplot(df7) +
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
		p <- ggplot(df7) +
			geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
			geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
			theme_bw() +
			theme(axis.text=element_text(color="black"),
				axis.ticks=element_line(color="black"))+
			scale_y_continuous(expand=c(0,0)) +
			scale_x_continuous("Ka/Ks") +
			scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
				"TE-like"="#E69F00","Unmethylated"="#CC79A7","Unclassified"="#999999"),
				breaks=c("Total","gbM","TE-like","Unmethylated","Unclassified"))
		ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
	}
}













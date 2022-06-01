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

#transposed genes
status=c("Aduranensis"="depleted",
	"Aipaensis"="depleted",
	"Alyrata"="depleted",
	"Athaliana"="depleted",
	"Atrichopoda"="depleted",
	"Bdistachyon"="depleted",
	"Boleracea"="depleted",
	"Brapa"="NS",
	"Bvulgaris"="depleted",
	"Cclementina"="depleted",
	"Cpapaya"="depleted",
	"Clanatus"="depleted",
	"Cmelo"="depleted",
	"Crubella"="depleted",
	"Csativus"="depleted",
	"Egrandis"="depleted",
	"Eguineensis"="NS",
	"Esalsugineum"="depleted",
	"Fvesca"="depleted",
	"Fxananassa"="enriched",
	"Gmax"="enriched",
	"Graimondii"="NS",
	"Ljaponicus"="depleted",
	"Macuminata"="enriched",
	"Mdomestica"="enriched",
	"Mesculenta"="NS",
	"Mguttatus"="depleted",
	"Mtruncatula"="depleted",
	"Osativa"="depleted",
	"Phallii"="depleted",
	"Ppersica"="depleted",
	"Ptrichocarpa"="enriched",
	"Pvirgatum"="NS",
	"Pvulgaris"="depleted",
	"Pxbretschneideri"="depleted",
	"Sbicolor"="depleted",
	"Sitalica"="depleted",
	"Slycopersicum"="enriched",
	"Stuberosum"="depleted",
	"Sviridis"="depleted",
	"Tcacao"="depleted",
	"Vvinifera"="enriched",
	"Zmays"="enriched")

Ks <- data.frame()
#Loop over each species
for(a in species){
	#Genes classified by mode of duplication
	path1 <-  paste(a,"/dupgen/results/kaks_results","/",a,".transposed.kaks",sep="")
	#Read in the files
	df1 <- read.table(path1,header=TRUE,sep="\t")[4]
	colnames(df1) <- c(a)
	df1 <- melt(df1)
	df1$median <- median(na.omit(df1$value))
	df1$status <- status[a]
	Ks <- rbind(Ks,df1)
}

write.csv(Ks,"../scripts/analyses_R/Ks.csv",quote=F)

p <- ggplot(Ks) + 
		geom_boxplot(aes(y=value,x=reorder(variable,mean),fill=status)) + 
		scale_fill_manual("TE-like genes",values=values2) + 
		scale_y_continuous("Ks",expand=c(0,0),limits=c(0,10)) + 
		xlab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))

ggsave("../../figures_tables/transposed_ks.pdf",p,device="pdf")

#tandem genes
status=c("Aduranensis"="",
	"Aipaensis"="",
	"Alyrata"="",
	"Athaliana"="",
	"Atrichopoda"="depleted",
	"Bdistachyon"="enriched",
	"Boleracea"="",
	"Brapa"="",
	"Bvulgaris"="",
	"Cclementina"="enriched",
	"Cpapaya"="",
	"Clanatus"="",
	"Cmelo"="",
	"Crubella"="",
	"Csativus"="",
	"Egrandis"="enriched",
	"Eguineensis"="enriched",
	"Esalsugineum"="",
	"Fvesca"="",
	"Fxananassa"="",
	"Gmax"="",
	"Graimondii"="",
	"Ljaponicus"="",
	"Macuminata"="NS",
	"Mdomestica"="",
	"Mesculenta"="",
	"Mguttatus"="",
	"Mtruncatula"="",
	"Osativa"="depleted",
	"Phallii"="depleted",
	"Ppersica"="",
	"Ptrichocarpa"="",
	"Pvirgatum"="depleted",
	"Pvulgaris"="",
	"Pxbretschneideri"="",
	"Sbicolor"="enriched",
	"Sitalica"="NS",
	"Slycopersicum"="",
	"Stuberosum"="",
	"Sviridis"="depleted",
	"Tcacao"="depleted",
	"Vvinifera"="NS",
	"Zmays"="depleted")

Ks <- data.frame()
#Loop over each species
for(a in species){
	#Genes classified by mode of duplication
	path1 <-  paste(a,"/dupgen/results/kaks_results","/",a,".tandem.kaks",sep="")
	#Read in the files
	df1 <- read.table(path1,header=TRUE,sep="\t")[4]
	colnames(df1) <- c(a)
	df1 <- melt(df1)
	df1$median <- median(na.omit(df1$value))
	df1$status <- status[a]
	Ks <- rbind(Ks,df1)
}

write.csv(Ks,"../scripts/analyses_R/Ks.csv",quote=F)

p <- ggplot(Ks) + 
		geom_boxplot(aes(y=value,x=reorder(variable,mean),fill=status)) + 
		scale_fill_manual("TE-like genes",values=values2) + 
		scale_y_continuous("Ks",expand=c(0,0),limits=c(0,10)) + 
		xlab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))

ggsave("../../figures_tables/tandem_ks.pdf",p,device="pdf")

#proximal genes
status=c("Aduranensis"="",
	"Aipaensis"="",
	"Alyrata"="",
	"Athaliana"="",
	"Atrichopoda"="depleted",
	"Bdistachyon"="",
	"Boleracea"="",
	"Brapa"="",
	"Bvulgaris"="",
	"Cclementina"="",
	"Cpapaya"="",
	"Clanatus"="",
	"Cmelo"="",
	"Crubella"="",
	"Csativus"="",
	"Egrandis"="",
	"Eguineensis"="",
	"Esalsugineum"="",
	"Fvesca"="",
	"Fxananassa"="",
	"Gmax"="",
	"Graimondii"="",
	"Ljaponicus"="",
	"Macuminata"="",
	"Mdomestica"="",
	"Mesculenta"="",
	"Mguttatus"="",
	"Mtruncatula"="",
	"Osativa"="",
	"Phallii"="",
	"Ppersica"="",
	"Ptrichocarpa"="",
	"Pvirgatum"="",
	"Pvulgaris"="",
	"Pxbretschneideri"="",
	"Sbicolor"="",
	"Sitalica"="",
	"Slycopersicum"="",
	"Stuberosum"="",
	"Sviridis"="",
	"Tcacao"="",
	"Vvinifera"="",
	"Zmays"="")

Ks <- data.frame()
#Loop over each species
for(a in species){
	#Genes classified by mode of duplication
	path1 <-  paste(a,"/dupgen/results/kaks_results","/",a,".proximal.kaks",sep="")
	#Read in the files
	df1 <- read.table(path1,header=TRUE,sep="\t")[4]
	colnames(df1) <- c(a)
	df1 <- melt(df1)
	df1$median <- median(na.omit(df1$value))
	df1$status <- status[a]
	Ks <- rbind(Ks,df1)
}

write.csv(Ks,"../scripts/analyses_R/Ks.csv",quote=F)

p <- ggplot(Ks) + 
		geom_boxplot(aes(y=value,x=reorder(variable,mean),fill=status)) + 
		scale_fill_manual("TE-like genes",values=values2) + 
		scale_y_continuous("Ks",expand=c(0,0),limits=c(0,10)) + 
		xlab("") + 
		theme_bw() + 
		theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))

ggsave("../../figures_tables/proximal_ks.pdf",p,device="pdf")

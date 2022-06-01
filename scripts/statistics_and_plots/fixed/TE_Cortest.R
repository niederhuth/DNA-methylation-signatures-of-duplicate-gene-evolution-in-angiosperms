library(ggplot2)
library(reshape2)
library(scales)
library(Hmisc)
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)


#This function provides a simple formatting of a correlation matrix
#into a table with 4 columns containing :
# Column 1 : row names (variable 1 for the correlation test)
# Column 2 : column names (variable 2 for the correlation test)
# Column 3 : the correlation coefficients
# Column 4 : the p-values of the correlations
flat_cor_mat <- function(cor_r, cor_p){
	cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
	cor_r <- gather(cor_r, column, cor, -1)
	cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
	cor_p <- gather(cor_p, column, p, -1)
	cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
	cor_p_matrix
}

#species = "Athaliana"

species = c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
	"Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
	"Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
	"Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
	"Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
	"Vvinifera","Zmays")

all_cor2 <- all_cor <- data.frame()

#### All-species TE-corr-loop ####
for( a in species){
	saving <- paste("../figures_tables/",a,sep="")
	if(!file.exists(saving)){
		dir.create(saving)
	}
	TE_path <- paste(a,"/methylpy/results/",a,"_TE_gene_distribution.tsv",sep="")
	TE_path2 <- paste(a,"/methylpy/results/",a,"_TE_gene_distribution_duplicates.tsv",sep="")

	#All Genes
	df1 <- read.table(TE_path, header = TRUE)
	df2 <- subset(df1,select = -c(Chr,Window))		# Removing columns Chr and Window
													# not needed for the correlation matrix
	cor_3 <- rcorr(as.matrix(df2), type = "pearson")  #Pearson correlations
	my_cor_matrix <- flat_cor_mat(cor_3$r, cor_3$P)
	my_cor_matrix$sps <- a

	write.csv(my_cor_matrix, file =(paste(saving,"/",a,"_TE_PearsonCorr.csv",sep="")))
	all_cor <- rbind(all_cor,as.data.frame(my_cor_matrix))

	df3 <- df2[,2:6]

	cor_4 <- rcorr(as.matrix(df3))  # unclassified is removed
	my_cor_matrix <- flat_cor_mat(cor_4$r, cor_4$P)

	# Making correlation plots for individual species for supplemental figure.
	p_mat <- cor_4$P
	test <- cor_4$r

	colnames(test) <- c("TEs","TE-nucleotides","gbM genes","teM genes","UnM genes")
	rownames(test) <- c("TEs","TE-nucleotides","gbM genes","teM genes","UnM genes")

	pdf(file = paste(saving,"/",a,"_TE_correlations_6x6.pdf", sep=""), width = 6, height = 6)
	col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
	corrplot(test, method = "color", col = col(20), type = "upper", tl.col = "darkblue", 
		tl.srt = 45, p.mat = p_mat, sig.level = 0.001)
	dev.off()

	#Duplicate Genes
	df1 <- read.table(TE_path2, header = TRUE)
	df2 <- subset(df1,select = -c(Chr,Window))		# Removing columns Chr and Window
													# not needed for the correlation matrix
	cor_3 <- rcorr(as.matrix(df2), type = "pearson")  #Pearson correlations
	my_cor_matrix2 <- flat_cor_mat(cor_3$r, cor_3$P)
	my_cor_matrix2$sps <- a

	write.csv(my_cor_matrix2,file=(paste(saving,"/",a,"_TE_PearsonCorr_duplicates.csv",sep="")))
	all_cor2 <- rbind(all_cor2,as.data.frame(my_cor_matrix2))

	df3 <- df2[,2:6]

	cor_4 <- rcorr(as.matrix(df3))  # unclassified is removed
	my_cor_matrix2 <- flat_cor_mat(cor_4$r, cor_4$P)

	# Making correlation plots for individual species for supplemental figure.
	p_mat <- cor_4$P
	test <- cor_4$r

	colnames(test) <- c("TEs","TE-nucleotides","gbM genes","teM genes","UnM genes")
	rownames(test) <- c("TEs","TE-nucleotides","gbM genes","teM genes","UnM genes")

	pdf(file = paste(saving,"/",a,"_TE_correlations_duplicates_6x6.pdf",sep=""),width=6,height=6)
	col <- colorRampPalette(c("#BB4444","#EE9988","#FFFFFF","#77AADD","#4477AA"))
	corrplot(test,method="color",col=col(20),type="upper",tl.col="darkblue", 
		tl.srt=45,p.mat=p_mat,sig.level=0.001)
	dev.off()
}
all_cor$p.adjust <- p.adjust(all_cor$p,method="BH")
all_cor2$p.adjust <- p.adjust(all_cor2$p,method="BH")

saving <- "../figures_tables/TE_correlations/"
if(!file.exists(saving)){
	dir.create(saving)
}
write.csv(all_cor,paste(saving,"All_TE_PearsonCorr.csv",sep=""),row.names=FALSE,quote=FALSE)
write.csv(all_cor2,paste(saving,"All_TE_PearsonCorr_duplicates.csv",sep=""),row.names=FALSE,
	quote=FALSE)



   
        
        



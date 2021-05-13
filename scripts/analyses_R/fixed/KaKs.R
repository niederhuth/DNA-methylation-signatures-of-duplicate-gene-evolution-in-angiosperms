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

setcolors <- c("Total"="black","gbM-gbM"="#F94144",
               "gbM-teM"="#F8961E","gbM-unM"="#F9C74F", "unM-unM" = "#90BE60", "unM-teM"="#43AA8B","teM-teM"="#277DA1")

setcolors2 <- c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1")

#Create empty dataframe for K-S test results
df13 <- df12 <- KaKs <- Ks <- data.frame()

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
	path4 <- paste("../figures_tables/", a,"/",a, "_duplicate_pair_met.csv",sep="")
	df4 <- read.csv(path4,header=TRUE)
	df9 <- data.frame()
	#Loop through each category of gene duplication
	for(b in c("wgd","proximal","dispersed","tandem","transposed")){
		#Read in appropriate kaks results
		path5a <- paste(a,"/dupgen/results/kaks_results","/",a,".",b,".kaks",sep="")
		path5u <- paste(a,"/dupgen/results-unique/kaks_results","/",a,".",b,".kaks",sep="")
		df5a <- read.table(path5a,header=T,sep="\t")
		tmp <- df5a[c(2,1,3,4,5,6)]
		colnames(tmp) <- colnames(df5a)
		df5a <- rbind(df5a,tmp)
		df5u <- read.table(path5u,header=T,sep="\t")
		tmp <- df5u[c(2,1,3,4,5,6)]
		colnames(tmp) <- colnames(df5u)
		df5u <- rbind(df5u,tmp)
		df6 <- merge(df4[df4$Duplication==b,],df5u,by=c("Duplicate.1","Duplicate.2"))
		df7 <- setdiff(df4[df4$Duplication==b,],df6[c(1:7)])
		df8 <- rbind(df6,merge(df7,df5a,by=c("Duplicate.1","Duplicate.2")))
		df6 <- merge(df4[df4$Duplication==b,],df5u,by=c("Duplicate.1","Duplicate.2"))
		df7 <- setdiff(df4[df4$Duplication==b,],df6[c(1:7)])
		df8 <- rbind(df6,merge(df7,df5a,by=c("Duplicate.1","Duplicate.2")))
	df9 <- rbind(df9,df8)
	}
	write.csv(df9,paste(path1,"/",a,"_kaks_values.csv",sep=""),row.names=FALSE,quote=FALSE)
	#Perform K-S test on each distribution
	for(c in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(nrow(df9[df9$Classification == c,]) != 0){
			tmp <- ks.test(df9[df9$Classification == c,]$Ks,
					sample(df9$Ks,size=nrow(df9[df9$Classification == c,]),replace=FALSE))
			Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
										D.statistic=tmp$stat,p.value=tmp$p.value))
			tmp <- ks.test(df9[df9$Classification == c,]$Ka.Ks,
					sample(df9$Ka.Ks,size=nrow(df9[df9$Classification == c,]),replace=FALSE))
			KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
										D.statistic=tmp$stat,p.value=tmp$p.value))
		} else {
			Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
										D.statistic=NA,p.value=NA))
			KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
										D.statistic=NA,p.value=NA))
		}
	}

	#remove duplicates with undetermined similarity
	df10 <- df9[df9$Similarity != "Undetermined",]
	#Add an order line for plots
	df10$order <- NA
	df10$order <- ifelse(df10$Classification=="gbM-gbM",1,df10$order)
	df10$order <- ifelse(df10$Classification=="gbM-teM",2,df10$order)
	df10$order <- ifelse(df10$Classification=="gbM-unM",3,df10$order)
	df10$order <- ifelse(df10$Classification=="unM-unM",4,df10$order)
	df10$order <- ifelse(df10$Classification=="unM-teM",5,df10$order)
	df10$order <- ifelse(df10$Classification=="teM-teM",6,df10$order)

	#Order the medians
	df10$Classification <- factor(df10$Classification, levels = c("gbM-gbM","gbM-teM","gbM-unM", "unM-unM","unM-teM","teM-teM"))
	df11 <- data.frame()
	for(b in c("gbM-gbM","gbM-teM","gbM-unM", "unM-unM","unM-teM","teM-teM")){
		df11 <- rbind(df11,data.frame(dup="SGD",Classification=b,median=median(na.omit(df10[df10$Classification==b & df10$Duplication != "wgd",]$Ks))))
	}
	df12 <- rbind(df12, data.frame(species=a,Duplication="SGD",data.frame(t(df11[order(df11$median),]$Classification))))
	df11 <- data.frame()
	for(b in c("gbM-gbM","gbM-teM","gbM-unM", "unM-unM","unM-teM","teM-teM")){
		df11 <- rbind(df11,data.frame(dup="WGD",Classification=b,median=median(na.omit(df10[df10$Classification==b & df10$Duplication == "wgd",]$Ks))))
	}
	df12 <- rbind(df12, data.frame(species=a,Duplication="WGD",data.frame(t(df11[order(df11$median),]$Classification))))
	df11 <- data.frame()
	for(b in c("gbM-gbM","gbM-teM","gbM-unM", "unM-unM","unM-teM","teM-teM")){
		df11 <- rbind(df11,data.frame(dup="SGD",Classification=b,median=median(na.omit(df10[df10$Classification==b & df10$Duplication != "wgd",]$Ka.Ks))))
	}
	df13 <- rbind(df13, data.frame(species=a,Duplication="SGD",data.frame(t(df11[order(df11$median),]$Classification))))
	df11 <- data.frame()
	for(b in c("gbM-gbM","gbM-teM","gbM-unM", "unM-unM","unM-teM","teM-teM")){
		df11 <- rbind(df11,data.frame(dup="WGD",Classification=b,median=median(na.omit(df10[df10$Classification==b & df10$Duplication == "wgd",]$Ka.Ks))))
	}
	df13 <- rbind(df13, data.frame(species=a,Duplication="WGD",data.frame(t(df11[order(df11$median),]$Classification))))

	### Lets Make Plots ###
	p <- ggplot(df10[df10$Duplication != "wgd",]) +
		geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
		geom_density(aes(x=Ks,color='Total'),size=1) +
		theme_bw() +
		theme(legend.position = "none") + 
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
		scale_color_manual(values = setcolors)
	ggsave(paste(path1,"/",a,"_SGDs_Ks_density.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication != "wgd",]) + 
		geom_boxplot(aes(y=Ks,x=reorder(Classification,order),color=Classification),size=1) + 
		theme_bw()+
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
		theme(axis.title.x = element_blank(),legend.position="none")+
		scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) +
		scale_color_manual(values=setcolors) +
		scale_y_continuous(expand=c(0.05,0.05),limits=c(0,5))
	ggsave(paste(path1,"/",a,"_SGDs_Ks_boxplot.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication != "wgd",]) +
		geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
		geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
		theme_bw() +
		theme(legend.position = "none") + 
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		geom_vline(xintercept=1, linetype="longdash", color="grey55") +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
		scale_color_manual(values = setcolors) 
	ggsave(paste(path1,"/",a,"_SGDs_KaKs_density.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication != "wgd",]) + 
		geom_boxplot(aes(y=Ka.Ks,x=reorder(Classification,order),color=Classification),size=1) + 
		theme_bw()+
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
		theme(axis.title.x = element_blank(),legend.position="none")+
		geom_hline(yintercept=1, linetype="longdash", color="grey55") +
		scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) +
		scale_color_manual(values=setcolors) +
		scale_y_continuous(expand=c(0.05,0.05),limits=c(0,2))
	ggsave(paste(path1,"/",a,"_SGDs_KaKs_boxplot.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication == "wgd",]) +
		geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
		geom_density(aes(x=Ks,color='Total'),size=1) +
		theme_bw() +
		theme(legend.position = "none") + 
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
		scale_color_manual(values = setcolors)
	ggsave(paste(path1,"/",a,"_WGDs_Ks_density.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication == "wgd",]) + 
		geom_boxplot(aes(y=Ks,x=reorder(Classification,order),color=Classification),size=1) + 
		theme_bw()+
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
		theme(axis.title.x = element_blank(),legend.position="none")+
		scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) +
		scale_color_manual(values=setcolors) +
		scale_y_continuous(expand=c(0.05,0.05),limits=c(0,5))
	ggsave(paste(path1,"/",a,"_WGDs_Ks_boxplot.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication == "wgd",]) +
		geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
		geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
		theme_bw() +
		theme(legend.position = "none") + 
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		geom_vline(xintercept=1, linetype="longdash", color="grey55") +
		scale_y_continuous(expand=c(0,0)) +
		scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
		scale_color_manual(values = setcolors)
	ggsave(paste(path1,"/",a,"_WGDs_KaKs_density.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

	p <- ggplot(df10[df10$Duplication == "wgd",]) + 
		geom_boxplot(aes(y=Ka.Ks,x=reorder(Classification,order),color=Classification),size=1) + 
		theme_bw()+
		theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
		theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
		theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
		theme(axis.title.x = element_blank(),legend.position="none")+
		geom_hline(yintercept=1, linetype="longdash", color="grey55") +
		scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) +
		scale_color_manual(values=setcolors) +
		scale_y_continuous(expand=c(0.05,0.05),limits=c(0,2))
	ggsave(paste(path1,"/",a,"_WGDs_KaKs_boxplot.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)

}

#Path for figures for all species
path6 <- "../figures_tables/kaks/"
if(!file.exists(path6)){
	dir.create(path6)
}
#Adjust p.values for K-S tests
Ks$p.adjust <- p.adjust(Ks$p.value,method="BH")
KaKs$p.adjust <- p.adjust(KaKs$p.value,method="BH")
#Output K-S results
write.csv(Ks,paste(path6,"Ks-distribution-test.csv",sep=""),row.names=FALSE,quote=FALSE)
write.csv(KaKs,paste(path6,"KaKs-distribution-test.csv",sep=""),row.names=FALSE,quote=FALSE)

write.csv(df12,paste(path6,"ks-median.csv",ep=""),row.names=FALSE,quote=FALSE)
write.csv(df13,paste(path6,"kaks-median.csv",sep=""),row.names=FALSE,quote=FALSE)

#Plot out Ks & KaKs medians
#Ks
df14 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df12[df12$Duplication=="SGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df14 <- rbind(df14,tmp)
}

p <- ggplot(df14) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand = c(0,0),limit=c(0,40)) + 
		scale_x_discrete("Lower Median Ks - Higher Median Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave(paste(path6,"SGD_median_Ks_order.pdf",sep=""),p,device="pdf")

df14 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df12[df12$Duplication=="WGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df14 <- rbind(df14,tmp)
}

p <- ggplot(df14) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand=c(0,0),limit=c(0,40)) + 
		scale_x_discrete("Lower Median Ks - Higher Median Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave(paste(path6,"WGD_median_Ks_order.pdf",sep=""),p,device="pdf")

#KaKs 
df14 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df13[df13$Duplication=="SGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df14 <- rbind(df14,tmp)
}

p <- ggplot(df14) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand=c(0,0),limit=c(0,40)) + 
		scale_x_discrete("Lower Median Ka/Ks - Higher Median Ka/Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave(paste(path6,"SGD_median_KaKs_order.pdf",sep=""),p,device="pdf")

df14 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df13[df13$Duplication=="WGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df14 <- rbind(df14,tmp)
}

p <- ggplot(df14) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand=c(0,0),limit=c(0,40)) + 
		scale_x_discrete("Lower Median Ka/Ks - Higher Median Ka/Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave(paste(path6,"WGD_median_KaKs_order.pdf",sep=""),p,device="pdf")







	











library(ggplot2)
library(reshape2)
library(scales)

species <- read.csv('../misc/genomes.csv',header=T)

for( a in species[species$methylC == "yes",]$species){
  path1 <- paste("../figures_tables/",a,sep="")
  if(!file.exists(path1)){
    dir.create(path1)
  }
}






	for(i in c("TE-like","gbM","Unmethylated","Unclassified")){
		ks.test(df7[df7$Classification == i,]$Ks,
			sample(df7$Ks,size=nrow(df7[df7$Classification == i,]),replace=FALSE))
		ks.test(df7[df7$Classification == i,]$Ks,
			df7[df7$Classification != i,]$Ks)
		ks.test(df7[df7$Classification == i,]$Ks,
			df7$Ks)
	}
}

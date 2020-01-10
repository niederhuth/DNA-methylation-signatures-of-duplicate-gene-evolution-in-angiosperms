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







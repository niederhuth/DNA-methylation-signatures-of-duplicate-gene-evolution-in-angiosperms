df1 <- read.table("dupgen/results-unique/classified_genes.tsv",header=TRUE,sep="\t")
df2 <- read.table("methylpy/results/classified_genes.txt",header=TRUE,sep="\t")[1,23]
df3 <- merge(df1,df2,by.x="Feature",by.y="Feature")
df4 <- as.data.frame(table(df3[2:3]))
df4$Perc <- NA
df4$Perc2 <- NA

for(i in df4$Duplication){
  x = sum(df4[df4$Duplication==i,]$Freq)
  df4[df4$Duplication==i,]$Perc = df4[df4$Duplication==i,]$Freq/x
}

for(i in df4$Classification){
  x = sum(df4[df4$Classification==i,]$Freq)
  df4[df4$Classification==i,]$Perc2 = df4[df4$Classification==i,]$Freq/x
}

p <- ggplot(df4) + 
      geom_bar(aes(x=Duplication,y=Perc2,fill=Classification),position="dodge",stat="identity")











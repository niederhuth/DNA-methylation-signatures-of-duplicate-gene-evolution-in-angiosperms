txdb <- makeTxDbFromGFF("Pvirgatum/ref/annotations/Pvirgatum.gff",format="gff")
introns <- intronsByTranscript(txdb,use.names=TRUE)
introns.tx <- data.frame(tx = 1:length(introns), exons = sapply(introns, length))

df1 <- read.table("Pvirgatum/dupgen/results-unique/Pvirgatum.wgd.pairs-unique",sep="\t",header=TRUE)[c(1,3)]
df2 <- read.table("Pvirgatum/methylpy/results/Pvirgatum_classified_genes.tsv",header=TRUE,sep="\t")

df1 <- rbind(df1,data.frame(Duplicate.1=df1$Duplicate.2,Duplicate.2=df1$Duplicate.1))

df3 <- merge(df1,df2[c(1,23)],by.x="Duplicate.2",by.y="Feature")
df3 <- merge(df3,df2[c(1,23)],by.x="Duplicate.1",by.y="Feature")
colnames(df3) <- c("Duplicate.1","Duplicate.2","Classification.1","Classification.2")

df4 <- na.omit(df3[df3$Classification.1=="gbM" & df3$Classification.2=="Unmethylated",])
df4 <- rbind(df4,na.omit(df3[df3$Classification.1=="gbM" & df3$Classification.2=="TE-like",]))
df4 <- rbind(df4,na.omit(df3[df3$Classification.1=="Unmethylated" & df3$Classification.2=="TE-like",]))

write.table(df4,"Pvirgatum/dupgen/results-unique/wgd_switch.tsv",
	quote=FALSE,sep="\t",row.names=FALSE)

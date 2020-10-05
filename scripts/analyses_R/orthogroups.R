library(ggplot2)
library(reshape2)
library(scales)

df1 <- read.table("methylpy/results/Athaliana_classified_genes.tsv",sep="\t",header=TRUE)
df1$Classification <- as.character(df1$Classification)
df1$Classification <- ifelse(is.na(df1$Classification),"Missing",df1$Classification)
df2 <- read.table("ref/mcscanx/Athaliana_orthogroups.tsv",header=FALSE,sep="\t")
df3 <- merge(df2,df1[c(1,23)],by.x="V1",by.y="Feature")
df4 <- data.frame(table(df3[c(2,3)]))
df5 <- df4[df4$Classification=="gbM",c(1,3)]
df5 <- merge(df5,df4[df4$Classification=="TE-like",c(1,3)],by="V2")
df5 <- merge(df5,df4[df4$Classification=="Unmethylated",c(1,3)],by="V2")
df5 <- merge(df5,df4[df4$Classification=="Unclassified",c(1,3)],by="V2")
df5 <- merge(df5,df4[df4$Classification=="Missing",c(1,3)],by="V2")
colnames(df5) <- c("Orthogroup","gbM","TE-like","Unmethylated","Unclassified","Missing")
df5$Total <- rowSums(df5[c(2:6)])
df6 <- melt(df5)


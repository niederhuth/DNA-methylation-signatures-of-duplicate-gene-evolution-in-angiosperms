library(ggplot2)

df1 <- read.table("Orthogroups.GeneCount.tsv",header=T,sep="\t",row.names=1)
df3 <- df2 <- subset(df1,select = -Total )
for(i in colnames(df2)){
  df2[i] <- ifelse(df1[i]==0,0,1)
}
df2$Total <- rowSums(df2)

for(i in colnames(df3)){
  df3[i] <- ifelse(df1[i]==1,1,0)
}
df3$Total <- rowSums(df3)/58

nrow(df2[df2$Total >= 51,])

df4 <- df3[row.names(df3) %in% row.names(df2[df2$Total >= 51,]),]
df5 <- df1[row.names(df1) %in%  row.names(df4[df4$Total >= 0.7,]),]


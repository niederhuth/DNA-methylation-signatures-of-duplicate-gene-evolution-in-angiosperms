df1 <- read.table("Athaliana/dupgen/results-unique/collinear_data.tsv",header=TRUE,sep="\t")
df1$collinearity <- gsub("-.*","",df1$Name)
df1$collinearity2 <- gsub(".*-","",df1$Name)
df2 <- merge(df1[df1$collinearity2==1,],df1[df1$collinearity2==2,],by="collinearity")



df4 <- read.table("Athaliana/methylpy/results/TE_intersections.tsv",header=T,sep="\t")
df4$Gene <- gsub("^.*Name=","",df4$Gene)
#colnames(df4[,1]) <- c("Feature")

df7 <- merge(merge(df6,df4,by.x="Duplicate.1",by.y="Gene"),df4,by.x="Duplicate.2",by.y="Gene")

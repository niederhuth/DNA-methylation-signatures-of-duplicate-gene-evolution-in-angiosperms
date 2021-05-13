df <- read.table("methylpy/results/Athaliana_classified_genes.tsv",
	header=T,sep="\t")[c("Feature","Classification")]
df2 <- read.table("dupgen/results-unique/classified_genes.tsv",header=T,sep="\t")
df3 <- merge(df,df2,by="Feature")

df4 <- read.csv("methylpy/results/all_Athaliana_var.csv",header=TRUE,
	stringsAsFactors=FALSE)
df4$Feature <- gsub(".Araport11.447","",df4$Feature)
df4[is.na(df4)] <- "Missing"
df5 <- data.frame(Feature=character(),gbM=numeric(),'TE-like'=numeric(),
	Unmethylated=numeric(),Unclassified=numeric(),Missing=numeric())
for(i in c(1:nrow(df))){
	df[i,]
	x <- as.data.frame(table(t(df4[i,])))
	df5 <- rbind(df5,data.frame(Feature=df4[i,]$Feature,
		gbM=ifelse(x[x[1]=="gbM",]$Freq,x[x[1]=="gbM",]$Freq,NA)[1],
		'TE-like'=ifelse(x[x[1]=="TE-like",]$Freq,x[x[1]=="TE-like",]$Freq,NA)[1],
		Unmethylated=ifelse(x[x[1]=="Unmethylated",]$Freq,
			x[x[1]=="Unmethylated",]$Freq,NA)[1],
		Unclassified=ifelse(x[x[1]=="Unclassified",]$Freq,
			x[x[1]=="Unclassified",]$Freq,NA)[1],
		Missing=ifelse(x[x[1]=="Missing",]$Freq,x[x[1]=="Missing",]$Freq,NA)[1]))
}
df5$Total <- rowSums(df5[2:5],na.rm=TRUE)
df5$p_gbM <- df5$gbM/df5$Total
df5$p_TE.like <- df5$TE.like/df5$Total
df5$p_Unmethylated <- df5$Unmethylated/df5$Total
df5$p_Unclassified <- df5$Unclassified/df5$Total

df6 <- merge(df3,df5,by="Feature")

df6$gbMgroup <- NA
df6$gbMgroup <- ifelse(df6$p_gbM > 0.75, 5, df6$gbMgroup)
df6$gbMgroup <- ifelse(df6$p_gbM < 0.75, 4, df6$gbMgroup)
df6$gbMgroup <- ifelse(df6$p_gbM < 0.5, 3, df6$gbMgroup)
df6$gbMgroup <- ifelse(df6$p_gbM < 0.25, 2, df6$gbMgroup)
df6$gbMgroup <- ifelse(is.na(df6$p_gbM), 1, df6$gbMgroup)

df6$TEgroup <- NA
df6$TEgroup <- ifelse(df6$p_TE.like > 0.75, 5, df6$TEgroup)
df6$TEgroup <- ifelse(df6$p_TE.like < 0.75, 4, df6$TEgroup)
df6$TEgroup <- ifelse(df6$p_TE.like < 0.5, 3, df6$TEgroup)
df6$TEgroup <- ifelse(df6$p_TE.like < 0.25, 2, df6$TEgroup)
df6$TEgroup <- ifelse(is.na(df6$p_TE.like), 1, df6$TEgroup)

df6$UnMgroup <- NA
df6$UnMgroup <- ifelse(df6$p_Unmethylated > 0.75, 5, df6$UnMgroup)
df6$UnMgroup <- ifelse(df6$p_Unmethylated < 0.75, 4, df6$UnMgroup)
df6$UnMgroup <- ifelse(df6$p_Unmethylated < 0.5, 3, df6$UnMgroup)
df6$UnMgroup <- ifelse(df6$p_Unmethylated < 0.25, 2, df6$UnMgroup)
df6$UnMgroup <- ifelse(is.na(df6$p_Unmethylated), 1, df6$UnMgroup)

wgd <- read.table("dupgen/results/kaks_results/Athaliana.wgd.kaks",header=T,sep="\t")
wgd$Duplication <- c("wgd")
tandem <- read.table("dupgen/results/kaks_results/Athaliana.tandem.kaks",header=T,sep="\t")
tandem$Duplication <- c("tandem")
proximal <- read.table("dupgen/results/kaks_results/Athaliana.proximal.kaks",header=T,sep="\t")
proximal$Duplication <- c("proximal")
transposed <- read.table("dupgen/results/kaks_results/Athaliana.transposed.kaks",header=T,sep="\t")
transposed$Duplication <- c("transposed")
dispersed <- read.table("dupgen/results/kaks_results/Athaliana.dispersed.kaks",header=T,sep="\t")
dispersed$Duplication <- c("dispersed")

df7 <- rbind(wgd,tandem,proximal,transposed,dispersed)
colnames(df7) <- c("Feature","Duplicate.2","Ka","Ks","Ka.Ks","P.Value","Duplication")
df8 <- df7[c(2,1,3,4,5,6,7)]
colnames(df8) <- colnames(df7)
df7 <- rbind(df7,df8)
rm(df8)

df8 <- data.frame()
y=0
for(i in df6[df6$Duplication != "unclassified" & df6$Duplication != "singletons",]$Feature){
	x <- df6[df6$Feature == i,]$Duplication
	if(nrow(df7[df7$Feature==i & df7$Duplication == x,]) != 0){
		df8 <- rbind(df8,df7[row.names(df7) == sample(row.names(df7[df7$Feature==i & df7$Duplication == x,]),1),c(1:6)])
	} else {
		df8 <- rbind(df8,data.frame(Feature=i,Duplicate.2=NA,Ka=NA,Ks=NA,Ka.Ks=NA,P.Value=NA))
		y = y + 1
	}
} 
print(paste("Missing: ",y,sep=""))
df9 <- data.frame()
y=0
for(i in df6[df6$Duplication != "unclassified" & df6$Duplication != "singletons",]$Feature){
	x <- df6[df6$Feature == i,]$Duplication
	if(nrow(df7[df7$Feature==i & df7$Duplication == x,]) != 0){
		df9 <- rbind(df9,df7[row.names(df7) == sample(row.names(df7[df7$Feature==i & df7$Ks == min(df7[df7$Feature==i & df7$Duplication == x,]$Ks),]),1),c(1:6)])
	} else {
		df9 <- rbind(df9,data.frame(Feature=i,Duplicate.2=NA,Ka=NA,Ks=NA,Ka.Ks=NA,P.Value=NA))
		y = y + 1
	}
} 
print(paste("Missing: ",y,sep=""))
df10 <- data.frame()
y=0
for(i in df6[df6$Duplication != "unclassified" & df6$Duplication != "singletons",]$Feature){
	x <- df6[df6$Feature == i,]$Duplication
	if(nrow(df7[df7$Feature==i & df7$Duplication == x,]) != 0){
		df10 <- rbind(df10,df7[row.names(df7) == sample(row.names(df7[df7$Feature==i & df7$Ks == max(df7[df7$Feature==i & df7$Duplication == x,]$Ks),]),1),c(1:6)])
	} else {
		df10 <- rbind(df10,data.frame(Feature=i,Duplicate.2=NA,Ka=NA,Ks=NA,Ka.Ks=NA,P.Value=NA))
		y = y + 1
	}
} 
print(paste("Missing: ",y,sep=""))
rm(y)

df11 <- merge(df6,df8,by="Feature")
labels <- c("0%","<25%","25-50%","50-75%",">75%")
ggplot(df11) + geom_density(aes(x=Ks,color=as.character(gbMgroup)),size=1) + 
	scale_color_discrete(labels=labels)
ggplot(df11) + geom_density(aes(x=Ks,color=as.character(TEgroup)),size=1) + 
	scale_color_discrete(labels=labels)
ggplot(df11) + geom_density(aes(x=Ks,color=as.character(UnMgroup)),size=1) + 
	scale_color_discrete(labels=labels)

ggplot(df11) + geom_boxplot(aes(y=Ks,x=as.character(gbMgroup)),size=1) + 
	scale_x_discrete(labels=labels)
ggplot(df11) + geom_boxplot(aes(y=Ks,x=as.character(TEgroup)),size=1) + 
	scale_x_discrete(labels=labels)
ggplot(df11) + geom_boxplot(aes(y=Ks,x=as.character(UnMgroup)),size=1) + 
	scale_x_discrete(labels=labels)

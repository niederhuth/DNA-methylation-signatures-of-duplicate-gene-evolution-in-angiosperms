library(ggplot2)


setcolors <- c("gbM-gbM"="#F94144","gbM-teM"="#F8961E","gbM-unM"="#F9C74F",
	"unM-unM"="#90BE60","unM-teM"="#43AA8B","teM-teM"="#277DA1")

#Ks 
df1 <- read.csv("ks-median.csv",header=TRUE,row.names=1)

df2 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df1[df1$Duplication=="SGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df2 <- rbind(df2,tmp)
}

p <- ggplot(df2) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand = c(0,0),limit=c(0,36)) + 
		scale_x_discrete("Lower Median Ks - Higher Median Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave("SGD_median_Ks_order.pdf")

df2 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df1[df1$Duplication=="WGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df2 <- rbind(df2,tmp)
}

p <- ggplot(df2) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand=c(0,0),limit=c(0,36)) + 
		scale_x_discrete("Lower Median Ks - Higher Median Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave("WGD_median_Ks_order.pdf")

#KaKs 
df1 <- read.csv("kaks-median.csv",header=TRUE,row.names=1)

df2 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df1[df1$Duplication=="SGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df2 <- rbind(df2,tmp)
}

p <- ggplot(df2) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand=c(0,0),limit=c(0,36)) + 
		scale_x_discrete("Lower Median Ka/Ks - Higher Median Ka/Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave("SGD_median_KaKs_order.pdf")

df2 <- data.frame()
for(a in c("X1","X2","X3","X4","X5","X6")){
	tmp <- data.frame(Var2=a,table(df1[df1$Duplication=="WGD",a]))
	for(b in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
		if(!(b %in% tmp$Var1)){
			tmp <- rbind(tmp,data.frame(Var2=a,Var1=b,Freq=0))
		}
	}
	df2 <- rbind(df2,tmp)
}

p <- ggplot(df2) + geom_bar(aes(x=Var2,y=Freq,fill=Var1),stat="identity",position="dodge") +
		theme_bw() + scale_y_continuous("Number of Species",expand=c(0,0),limit=c(0,36)) + 
		scale_x_discrete("Lower Median Ka/Ks - Higher Median Ka/Ks",labels=c(1,2,3,4,5,6)) +
		scale_fill_manual("Methylation Pairs",values=setcolors)
ggsave("WGD_median_KaKs_order.pdf")

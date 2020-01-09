library(ggplot2)
library(reshape2)
library(scales)

species <- read.csv('../misc/genomes.csv',header=T)

for( a in species[species$methylC == "yes",]$species){
  path1 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
  path2 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
  df1 <- read.table(path1,header=TRUE,sep="\t")
  df2 <- read.table(path2,header=TRUE,sep="\t")[1,23]
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
  
}

#Test
species="Athaliana"

for( a in species){
  path1 <- "../../figures_tables/"
  dir.create(paste(path1,a,sep=""))
  path2 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
  path3 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
  df1 <- read.table(path2,header=TRUE,sep="\t")
  df2 <- read.table(path3,header=TRUE,sep="\t")
  df3 <- merge(df1,df2[,c(1,23)],by.x="Feature",by.y="Feature")
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
    geom_bar(aes(x=Duplication,y=Perc2,fill=Classification),position="dodge",stat="identity") +
    theme_bw() +
    scale_y_continuous(labels=percent,expand=c(0,0),limits=c(0,0.35))
  ggsave(paste(path1,a,"/",a,"_test1.pdf",sep=""),p,device="pdf")
  
  df8 <- data.frame(Classification.x.y=c("gbM-gbM","gbM-TE-like","gbM-Unclassified","gbM-Unmethylated","TE-like-TE-like",
                                "TE-like-Unclassified","Unclassified-Unclassified","Unmethylated-TE-like",
                                "Unmethylated-Unclassified","Unmethylated-Unmethylated"))
  df9 <- data.frame(Change=c("Identical","Different"))
  
  for(b in c("wgd","tandem","proximal","dispersed")){
    path4 <- paste(a,"/dupgen/results-unique/",a,".",b,".pairs-unique",sep="")
    df5 <- read.table(path4,header=TRUE,sep="\t")[,c(1,3)]
    df6 <- merge(df5,df2[,c(1,6,11,16,23)],by.x="Duplicate.1",by.y="Feature")
    df7 <- na.omit(merge(df6,df2[,c(1,6,11,16,23)],by.x="Duplicate.2",by.y="Feature"))
    df7$Change <- ifelse(df7$Classification.x==df7$Classification.y,"Identical","Different")
    df7$Classification.x.y <- "NA"
    
    for(row in 1:nrow(df7)){
      if(df7[row,]$Classification.x == "NA" || df7[row,]$Classification.y == "NA"){
        df7[row,]$Classification.x.y = df7[row,]$Classification.x.y
      } else if(df7[row,]$Classification.x == df7[row,]$Classification.y){
        df7[row,]$Classification.x.y = paste(df7[row,]$Classification.x,df7[row,]$Classification.y,sep="-")
      } else if(df7[row,]$Classification.x == "gbM" & df7[row,]$Classification.y == "Unclassified" || df7[row,]$Classification.x == "Unclassified" & df7[row,]$Classification.y == "gbM" ){
        df7[row,]$Classification.x.y = "gbM-Unclassified"
      }  else if(df7[row,]$Classification.x == "gbM" & df7[row,]$Classification.y == "TE-like" || df7[row,]$Classification.x == "TE-like" & df7[row,]$Classification.y == "gbM" ){
        df7[row,]$Classification.x.y = "gbM-TE-like"
      }  else if(df7[row,]$Classification.x == "gbM" & df7[row,]$Classification.y == "Unmethylated" || df7[row,]$Classification.x == "Unmethylated" & df7[row,]$Classification.y == "gbM" ){
        df7[row,]$Classification.x.y = "gbM-Unmethylated"
      } else if(df7[row,]$Classification.x == "Unmethylated" & df7[row,]$Classification.y == "TE-like" || df7[row,]$Classification.x == "TE-like" & df7[row,]$Classification.y == "Unmethylated" ){
        df7[row,]$Classification.x.y = "Unmethylated-TE-like"
      } else if(df7[row,]$Classification.x == "Unmethylated" & df7[row,]$Classification.y == "Unclassified" || df7[row,]$Classification.x == "Unclassified" & df7[row,]$Classification.y == "Unmethylated" ){
        df7[row,]$Classification.x.y = "Unmethylated-Unclassified"
      } else if(df7[row,]$Classification.x == "TE-like" & df7[row,]$Classification.y == "Unclassified" || df7[row,]$Classification.x == "Unclassified" & df7[row,]$Classification.y == "TE-like" ){
        df7[row,]$Classification.x.y = "TE-like-Unclassified"
      } else {
        df7[row,]$Classification.x.y = df7[row,]$Classification.x.y
      }
    }
    
    df8[,b] <- data.frame(table(df7$Classification.x.y))$Freq
    df9[,b] <- data.frame(table(df7$Change))$Freq
    
    p <- ggplot(df7) +
      geom_point(aes(x=log2((CG_Weighted_mC.y*100)),y=log2((CG_Weighted_mC.x*100)),color=Classification.x.y)) +
      theme_bw() +
      xlab("Duplicate 1") +
      ylab("Duplicate 2") + 
      geom_smooth(aes(x=log2((CG_Weighted_mC.y*100)),y=log2((CG_Weighted_mC.x*100))),method=lm)
    ggsave(paste(path1,a,"/",a,"_",b,"_CG.pdf",sep=""),p,device="pdf")
    
    p <- ggplot(df7) +
      geom_point(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100)),color=Classification.x.y)) +
      theme_bw() +
      xlab("Duplicate 1") +
      ylab("Duplicate 2") + 
      geom_smooth(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100))),method=lm)
    ggsave(paste(path1,a,"/",a,"_",b,"_CHG.pdf",sep=""),p,device="pdf")
  
    p <- ggplot(df7) +
      geom_point(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100)),color=Classification.x.y)) +
      theme_bw() +
      xlab("Duplicate 1") +
      ylab("Duplicate 2") + 
      geom_smooth(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100))),method=lm)
    ggsave(paste(path1,a,"/",a,"_",b,"_CHH.pdf",sep=""),p,device="pdf")
  }    
  
  path4 <- paste(a,"/dupgen/results-unique/",a,".transposed.pairs-unique",sep="")
  b="transposed"
  df5 <- read.table(path4,header=TRUE,sep="\t")[,c(1,3)]
  df6 <- merge(df5,df2[,c(1,6,11,16,23)],by.x="Transposed",by.y="Feature")
  df7 <- merge(df6,df2[,c(1,6,11,16,23)],by.x="Parental",by.y="Feature")
  df7$Change <- ifelse(df7$Classification.x==df7$Classification.y,"Identical","Different")
  
  df7$Classification.x <- ifelse(is.na(df7$Classification.x),"NA",df7$Classification.x)
  df7$Classification.y <- ifelse(is.na(df7$Classification.y),"NA",df7$Classification.y)
  
  df7$Classification.x.y <- "NA"
  
  for(row in 1:nrow(df7)){
    if(df7[row,]$Classification.x == "NA" || df7[row,]$Classification.y == "NA"){
      df7[row,]$Classification.x.y = "NA"
    } else if(df7[row,]$Classification.x == df7[row,]$Classification.y){
      df7[row,]$Classification.x.y = paste(df7[row,]$Classification.x,df7[row,]$Classification.y,sep="-")
    } else if(df7[row,]$Classification.x == "gbM" & df7[row,]$Classification.y == "Unclassified" || df7[row,]$Classification.x == "Unclassified" & df7[row,]$Classification.y == "gbM" ){
      df7[row,]$Classification.x.y = "gbM-Unclassified"
    }  else if(df7[row,]$Classification.x == "gbM" & df7[row,]$Classification.y == "TE-like" || df7[row,]$Classification.x == "TE-like" & df7[row,]$Classification.y == "gbM" ){
      df7[row,]$Classification.x.y = "gbM-TE-like"
    }  else if(df7[row,]$Classification.x == "gbM" & df7[row,]$Classification.y == "Unmethylated" || df7[row,]$Classification.x == "Unmethylated" & df7[row,]$Classification.y == "gbM" ){
      df7[row,]$Classification.x.y = "gbM-Unmethylated"
    } else if(df7[row,]$Classification.x == "Unmethylated" & df7[row,]$Classification.y == "TE-like" || df7[row,]$Classification.x == "TE-like" & df7[row,]$Classification.y == "Unmethylated" ){
      df7[row,]$Classification.x.y = "Unmethylated-TE-like"
    } else if(df7[row,]$Classification.x == "Unmethylated" & df7[row,]$Classification.y == "Unclassified" || df7[row,]$Classification.x == "Unclassified" & df7[row,]$Classification.y == "Unmethylated" ){
      df7[row,]$Classification.x.y = "Unmethylated-Unclassified"
    } else if(df7[row,]$Classification.x == "TE-like" & df7[row,]$Classification.y == "Unclassified" || df7[row,]$Classification.x == "Unclassified" & df7[row,]$Classification.y == "TE-like" ){
      df7[row,]$Classification.x.y = "TE-like-Unclassified"
    } else {
      df7[row,]$Classification.x.y = df7[row,]$Classification.x.y
    }
  }
  
  df8[,b] <- data.frame(table(df7$Classification.x.y))$Freq
  df9[,b] <- data.frame(table(df7$Change))$Freq
  
  p <- ggplot(df7) +
    geom_point(aes(x=log2((CG_Weighted_mC.y*100)),y=log2((CG_Weighted_mC.x*100)),color=Classification.x.y)) +
    theme_bw() +
    xlab("Parental") +
    ylab("Transposed") +   
    geom_smooth(aes(x=log2((CG_Weighted_mC.y*100)),y=log2((CG_Weighted_mC.x*100))),method=lm)
  ggsave(paste(path1,a,"/",a,"_",b,"_CG.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df7) +
    geom_point(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100)),color=Classification.x.y)) +
    theme_bw() +
    xlab("Parental") +
    ylab("Transposed") + 
    geom_smooth(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100))),method=lm)
  ggsave(paste(path1,a,"/",a,"_",b,"_CHG.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df7) +
    geom_point(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100)),color=Classification.x.y)) +
    theme_bw() +
    xlab("Parental") +
    ylab("Transposed") + 
    geom_smooth(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100))),method=lm)
  ggsave(paste(path1,a,"/",a,"_",b,"_CHH.pdf",sep=""),p,device="pdf")
  
  
  
  df9 <- melt(df8)
  p <- ggplot(df9) +
        geom_bar(aes(x=Classification.x.y,y=value,fill=variable),position="dodge",stat="identity") +
        theme_bw()
}







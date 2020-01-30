library(ggplot2)
library(reshape2)
library(scales)

#Test
species="Athaliana"

for( a in species){
  path1 <- paste("../../figures_tables/",a,sep="")
  if(!file.exists(path1)){
    dir.create(path1)
  }
  path2 <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
  path3 <- paste(a,"/methylpy/results/",a,"_classified_genes.tsv",sep="")
  df1 <- read.table(path2,header=TRUE,sep="\t")
  df2 <- read.table(path3,header=TRUE,sep="\t")
  df3 <- merge(df1,df2[,c(1,23)],by.x="Feature",by.y="Feature")
  df4 <- as.data.frame(table(df3[2:3]))
  df4$Perc2 <- df4$Perc <- NA
  
  for(i in df4$Duplication){
    x = sum(df4[df4$Duplication==i,]$Freq)
    df4[df4$Duplication==i,]$Perc = df4[df4$Duplication==i,]$Freq/x
  }
  
  for(i in df4$Classification){
    x = sum(df4[df4$Classification==i,]$Freq)
    df4[df4$Classification==i,]$Perc2 = df4[df4$Classification==i,]$Freq/x
  }
 
  p <- ggplot(df4) + 
    geom_bar(aes(x=Duplication,y=Perc,fill=Classification),position="dodge",stat="identity") +
    theme_bw() +
    scale_fill_discrete("Methylation Classification") +
    scale_y_continuous(labels=percent,expand=c(0,0)) +
    xlab("Duplication Type") +
    ylab("Percentage of Genes")
  #ggsave(paste(path1,"/",a,"_test1.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df4) + 
    geom_bar(aes(x=Classification,y=Perc2,fill=Duplication),position="dodge",stat="identity") +
    theme_bw() +
    scale_fill_discrete("Duplication Type") +
    scale_y_continuous(labels=percent,expand=c(0,0)) +
    xlab("Methylation Classification") +
    ylab("Percentage of Genes")
  #ggsave(paste(path1,"/",a,"_test2.pdf",sep=""),p,device="pdf")
  
  for(d in unique(df4$Duplication)){
   p <- ggplot(df4[df4$Duplication==d,],aes(x="",y=Perc,fill=Classification)) +
          geom_bar(width = 1, stat = "identity")+
          coord_polar("y", start=0) +
          theme_void() +
          theme(axis.text.x=element_blank()) +
          geom_text(aes(y = Perc/4 + c(0, cumsum(Perc)[-length(Perc)]), 
                    label = percent(sort(Perc))), size=5)
    #ggsave(paste(path1,"/",a,"_",d,"_pie.pdf",sep=""),p,device="pdf")
  }

  df8 <- data.frame(Classification.x.y=c("gbM-gbM","gbM-TE-like","gbM-Unclassified","gbM-Unmethylated","TE-like-TE-like",
                                "TE-like-Unclassified","Unclassified-Unclassified","Unmethylated-TE-like",
                                "Unmethylated-Unclassified","Unmethylated-Unmethylated"))
  df9 <- data.frame(Change=c("Identical","Different"))
  
  df19 <- data.frame()
  for(b in c("wgd","proximal","dispersed","tandem")){
    path4 <- paste(a,"/dupgen/results-unique/",a,".",b,".pairs-unique",sep="")
    df5 <- read.table(path4,header=TRUE,sep="\t")[,c(1,3)]
    df6 <- merge(df5,df2[,c(1,6,11,16,23)],by.x="Duplicate.1",by.y="Feature")
    df7 <- na.omit(merge(df6,df2[,c(1,6,11,16,23)],by.x="Duplicate.2",by.y="Feature"))
    df7$Change <- ifelse(df7$Classification.x==df7$Classification.y,"Identical","Different")
    df7$Classification.x.y <- "NA"
    
    for(row in 1:nrow(df7)){
      if(df7[row,]$Classification.x == df7[row,]$Classification.y){
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
    #ggsave(paste(path1,"/",a,"_",b,"_CG.pdf",sep=""),p,device="pdf")
    
    p <- ggplot(df7) +
      geom_point(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100)),color=Classification.x.y)) +
      theme_bw() +
      xlab("Duplicate 1") +
      ylab("Duplicate 2") + 
      geom_smooth(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100))),method=lm)
    #ggsave(paste(path1,"/",a,"_",b,"_CHG.pdf",sep=""),p,device="pdf")
  
    p <- ggplot(df7) +
      geom_point(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100)),color=Classification.x.y)) +
      theme_bw() +
      xlab("Duplicate 1") +
      ylab("Duplicate 2") + 
      geom_smooth(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100))),method=lm)
    #ggsave(paste(path1,"/",a,"_",b,"_CHH.pdf",sep=""),p,device="pdf")
    
    #kaks
    path5 <- paste(a,"/dupgen/results-unique/kaks_results","/",a,".",b,".kaks",sep="")
    df10 <- read.table(path5,header=T,sep="\t")
    df11 <- merge(df7,df10[,c(1,3:5)],by.x="Duplicate.1",by.y="Duplicate.1")
    
    p <- ggplot(df11) +
      geom_boxplot(aes(x=Classification.x.y,y=Ka.Ks)) +
      theme_bw() +
      xlab("Methylation Classification") +
      ylab("Ka/Ks")
    #ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
    
    #GC content stuff
    path8 <- paste(a,"/ref/mcscanx/",a,"-gc123.tsv",sep="")
    df16 <- read.table(path8,header=T,sep="\t")
    df11 <- merge(df11,df16,by.x="Duplicate.1",by.y="Transcript")
    df11 <- merge(df11,df16,by.x="Duplicate.2",by.y="Transcript")
    
    tmp <- df11[,c(2:6,11,12,15:19)]
    tmp2 <- df11[,c(1,7:12,15,20:23)]
    colnames(tmp2) <- colnames(tmp) <- c("Gene","CG_Weighted_mC","CHG_Weighted_mC","CHH_Weighted_mC",
                                         "Classification","Change","Classification.x.y","Ka.Ks","GC",
                                         "GC1","GC2","GC3")
    df17 <- rbind(tmp,tmp2)
    
    p <- ggplot(df17[df17$Change=="Different",]) + 
      geom_boxplot(aes(x=Classification.x.y,fill=Classification,y=GC),position="dodge") 
    
    p <- ggplot(df17[df17$Change=="Different",]) + 
      geom_boxplot(aes(x=Classification.x.y,fill=Classification,y=GC3),position="dodge")
    
    #test if different
    for(c in unique(df11[df11$Change=="Different",]$Classification.x.y)){
      d <- gsub("-.*","",c)
      tmp <- df11[df11$Classification.x.y==c & df11$Classification.x==d,c(16:23)]
      tmp2 <- df11[df11$Classification.x.y==c & df11$Classification.x!=d,c(20:23,16:19)]
      colnames(tmp2) <- colnames(tmp)
      df18 <- rbind(tmp,tmp2)
      df19 <- rbind(df19,data.frame(Duplication_type=b,
                                    mC_Class_Change=c,
                                    GC=t.test(df18$GC.x,df18$GC.y,paired=T)$p.value,
                                    GC1=t.test(df18$GC1.x,df18$GC1.y,paired=T)$p.value,
                                    GC2=t.test(df18$GC2.x,df18$GC2.y,paired=T)$p.value,
                                    GC3=t.test(df18$GC3.x,df18$GC3.y,paired=T)$p.value))
    }
    
    #test over/under representation
    df20 <- data.frame()
    for(e in 1:10000){
      df20 <- rbind(df20,as.data.frame(table(df2[df2$Feature %in% sample(df2$Feature,nrow(df11)*2,
                                                            replace = F),"Classification"]))$Freq)
    }
    colnames(df20) <- c("gbM","TE-like","Unclassified","Unmethylated")
    print(b)
    for(f in c('gbM','TE-like','Unclassified','Unmethylated')){
      print(f)
      df21 <- data.frame(table(df11$Classification.x)+table(df11$Classification.y))
      g <- df21[df21$Var1 == f,]$Freq
      if(g < mean(df20[,f])){
        h <- pnorm(g,mean=mean(df20[,f]),sd=sd(df20[,f]),lower.tail=T)*2
        print(h)
      } else {
        h <- pnorm(g,mean=mean(df20[,f]),sd=sd(df20[,f]),lower.tail=F)*2
        print(h)
      }  
    }
  }    
  
  #Tandem Gene Orientation
  path7 <- paste(a,"/dupgen/results-unique/orientation.tsv",sep="")
  df15 <- read.table(path7,header=F,sep="\t")
  df11$Duplicate.1.strand <- df11$Duplicate.2.strand <- NA
  for(row in 1:nrow(df11)){
    df11[row,]$Duplicate.1.strand <- as.character(df15[df15$V1==as.character(df11[row,]$Duplicate.1),]$V2)
    df11[row,]$Duplicate.2.strand <- as.character(df15[df15$V1==as.character(df11[row,]$Duplicate.2),]$V2)
  }
  df11$strand.switch <- ifelse(df11$Duplicate.1.strand == df11$Duplicate.2.strand,"Identical","Opposite")
  df11$strand.switch2 <- ifelse(df11$Duplicate.1.strand == "+" & df11$Duplicate.2.strand == "-",
                                "Inverted","Not-Inverted")
  
  #Transposed Genes
  path4 <- paste(a,"/dupgen/results-unique/",a,".transposed.pairs-unique",sep="")
  b="transposed"
  df5 <- read.table(path4,header=TRUE,sep="\t")[,c(1,3)]
  df6 <- merge(df5,df2[,c(1,6,11,16,23)],by.x="Transposed",by.y="Feature")
  df7 <- na.omit(merge(df6,df2[,c(1,6,11,16,23)],by.x="Parental",by.y="Feature"))
  df7$Change <- ifelse(df7$Classification.x==df7$Classification.y,"Identical","Different")
  df7$Classification.x.y <- "NA"
  
  for(row in 1:nrow(df7)){
    if(df7[row,]$Classification.x == df7[row,]$Classification.y){
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
  #ggsave(paste(path1,"/",a,"_",b,"_CG.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df7) +
    geom_point(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100)),color=Classification.x.y)) +
    theme_bw() +
    xlab("Parental") +
    ylab("Transposed") + 
    geom_smooth(aes(x=log2((CHG_Weighted_mC.y*100)),y=log2((CHG_Weighted_mC.x*100))),method=lm)
  #ggsave(paste(path1,"/",a,"_",b,"_CHG.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df7) +
    geom_point(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100)),color=Classification.x.y)) +
    theme_bw() +
    xlab("Parental") +
    ylab("Transposed") + 
    geom_smooth(aes(x=log2((CHH_Weighted_mC.y*100)),y=log2((CHH_Weighted_mC.x*100))),method=lm)
  #ggsave(paste(path1,"/",a,"_",b,"_CHH.pdf",sep=""),p,device="pdf")
  
  df8 <- melt(df8)
  df9 <- melt(df9)
  
  df8$Perc2 <- df8$Perc <- NA
  df9$Perc2 <- df9$Perc <- NA

  for(i in df8$Classification.x.y){
    x = sum(df8[df8$Classification.x.y==i,]$value)
    df8[df8$Classification.x.y==i,]$Perc = df8[df8$Classification.x.y==i,]$value/x
  }
  
  for(i in unique(df8$variable)){
    x = sum(df8[df8$variable==i,]$value)
    df8[df8$variable==i,]$Perc2 = df8[df8$variable==i,]$value/x
  }
  
  for(i in df9$Change){
    x = sum(df9[df9$Change==i,]$value)
    df9[df9$Change==i,]$Perc = df9[df9$Change==i,]$value/x
  }
  
  for(i in unique(df9$variable)){
    x = sum(df9[df9$variable==i,]$value)
    df9[df9$variable==i,]$Perc2 = df9[df9$variable==i,]$value/x
  }
  
  p <- ggplot(df8) +
    geom_bar(aes(x=Classification.x.y,y=Perc,fill=variable),position="dodge",stat="identity") +
    theme_bw() +
    scale_y_continuous(labels=percent,expand=c(0,0)) +
    scale_fill_discrete("Duplication Type") +
    xlab("Methylation Classification") +
    ylab("Percentage of Genes")
  #ggsave(paste(path1,"/",a,"_gene_pairs_class1.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df8) +
    geom_bar(aes(x=variable,y=Perc2,fill=Classification.x.y),position="dodge",stat="identity") +
    theme_bw() +
    scale_fill_discrete("Methylation Classification") +
    scale_y_continuous(labels=percent,expand=c(0,0)) +
    xlab("Duplication Type") +
    ylab("Percentage of Genes")
  #ggsave(paste(path1,"/",a,"_gene_pairs_class2.pdf",sep=""),p,device="pdf")

  p <- ggplot(df9) +
    geom_bar(aes(x=Change,y=Perc,fill=variable),position="dodge",stat="identity") +
    theme_bw() +
    scale_fill_discrete("Duplication Type") +
    scale_y_continuous(labels=percent,expand=c(0,0)) +
    xlab("Methylation Classification") +
    ylab("Percentage of Genes")
  #ggsave(paste(path1,"/",a,"_gene_pairs_change1.pdf",sep=""),p,device="pdf")
  
  p <- ggplot(df9) +
    geom_bar(aes(x=variable,y=Perc2,fill=Change),position="dodge",stat="identity") +
    theme_bw() +
    scale_fill_discrete("Methylation Classification") +
    scale_y_continuous(labels=percent,expand=c(0,0)) +
    xlab("Duplication Type") +
    ylab("Percentage of Genes")
  #ggsave(paste(path1,"/",a,"_gene_pairs_change2.pdf",sep=""),p,device="pdf")
  
  #kaks
  path5 <- paste(a,"/dupgen/results-unique/kaks_results","/",a,".",b,".kaks",sep="")
  df10 <- read.table(path5,header=T,sep="\t")
  df11 <- merge(df7,df10[,c(1,3:5)],by.x="Transposed",by.y="Duplicate.1")
  
  p <- ggplot(df11) +
    geom_boxplot(aes(x=Classification.x.y,y=Ka.Ks)) +
    theme_bw() +
    xlab("Methylation Classification") +
    ylab("Ka/Ks")
  #ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
  
  #GC stuff
  path8 <- paste(a,"/ref/mcscanx/",a,"-gc123.tsv",sep="")
  df16 <- read.table(path8,header=T,sep="\t")
  df11 <- merge(df11,df16,by.x="Transposed",by.y="Transcript")
  df11 <- merge(df11,df16,by.x="Parental",by.y="Transcript")
  
  tmp <- df11[,c(2:6,11,12,15:19)]
  tmp$Copy <- "Transposed"
  tmp2 <- df11[,c(1,7:12,15,20:23)]
  tmp2$Copy <- "Parental"
  colnames(tmp2) <- colnames(tmp) <- c("Gene","CG_Weighted_mC","CHG_Weighted_mC","CHH_Weighted_mC",
                                       "Classification","Change","Classification.x.y","Ka.Ks","GC",
                                       "GC1","GC2","GC3","Copy")
  df17 <- rbind(tmp,tmp2)
  
  p <- ggplot(df17[df17$Change=="Different",]) + 
    geom_boxplot(aes(x=Classification.x.y,fill=Classification,y=GC),position="dodge") 
  
  p <- ggplot(df17[df17$Change=="Different",]) + 
    geom_boxplot(aes(x=Classification.x.y,fill=Classification,y=GC3),position="dodge")
  
  #test if different
  for(c in unique(df11[df11$Change=="Different",]$Classification.x.y)){
    d <- gsub("-.*","",c)
    tmp <- df11[df11$Classification.x.y==c & df11$Classification.x==d,c(16:23)]
    tmp2 <- df11[df11$Classification.x.y==c & df11$Classification.x!=d,c(20:23,16:19)]
    colnames(tmp2) <- colnames(tmp)
    df18 <- rbind(tmp,tmp2)
    df19 <- rbind(df19,data.frame(Duplication_type=b,
                                  mC_Class_Change=c,
                                  GC=t.test(df18$GC.x,df18$GC.y,paired=T)$p.value,
                                  GC1=t.test(df18$GC1.x,df18$GC1.y,paired=T)$p.value,
                                  GC2=t.test(df18$GC2.x,df18$GC2.y,paired=T)$p.value,
                                  GC3=t.test(df18$GC3.x,df18$GC3.y,paired=T)$p.value))
  }

  #Transposed Copy Age
  path6 <- paste(a,"/mcscanx/results/Athaliana.transposed_epoch.pairs",sep="")
  if(file.exists(path6)){
    df12 <- read.table(path6,header=T,sep="\t")[,c(1,6,7)]
    df13 <- merge(df11,df12,by.x="Transposed",by.y="Transposed")
    p <- ggplot(df13) + 
      geom_boxplot(aes(x=reorder(epoch_species,epoch),y=Ka.Ks,fill=Classification.x),position="dodge") +
      theme_bw() +
      scale_fill_discrete("Methylation Classification") +
      xlab("Since Divergence From") +
      ylab("Ka/Ks")
    #ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
    
    p <- ggplot(df13) + 
      geom_bar(aes(x=reorder(epoch_species,epoch),fill=Classification.x),position="dodge") +
      theme_bw() +
      scale_fill_discrete("Methylation Classification") +
      xlab("Since Divergence From") +
      ylab("Number of Genes")
    #ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
    
    df14 <- data.frame()
    for(i in unique(df13$Classification.x)){
      df14 <- rbind(df14,data.frame(Classification=i,table(df13[df13$Classification.x==i,]$epoch_species)))
    }
    df14$Perc <- NA
    for(i in unique(df14$Classification)){
      x=sum(df14[df14$Classification==i,]$Freq)
      df14[df14$Classification==i,]$Perc = df14[df14$Classification==i,]$Freq/x
    }
    df14 <- merge(df14,unique(df13[,c("epoch","epoch_species")]),by.x="Var1",by.y="epoch_species")
    df14$Perc2 <- NA
    for(i in unique(df14$Var1)){
      x=sum(df14[df14$Var1==i,]$Freq)
      df14[df14$Var1==i,]$Perc2 = df14[df14$Var1==i,]$Freq/x
    }
    
    p <- ggplot(df14) + 
      geom_bar(aes(x=reorder(Var1,epoch),fill=Classification,y=Perc2),position="dodge",stat="identity") +
      theme_bw() +
      scale_fill_discrete("Methylation Classification") +
      xlab("Since Divergence From") +
      ylab("Percentage of Genes")
    #ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")
    
    p <- ggplot(df14) + 
      geom_bar(aes(x=Classification,fill=reorder(Var1,epoch),y=Perc),position="dodge",stat="identity") +
      theme_bw() +
      scale_fill_discrete("Since Divergence From") +
      xlab("Methylation Classification") +
      ylab("Percentage of Genes")
    #ggsave(paste(path1,"/",a,"_",b,"_KaKs.pdf",sep=""),p,device="pdf")    
  } 
  
  #test over/under representation
  df20 <- data.frame()
  for(e in 1:10000){
    df20 <- rbind(df20,as.data.frame(table(df2[df2$Feature %in% sample(df2$Feature,nrow(df11)*2,
                                                                       replace = F),"Classification"]))$Freq)
  }
  colnames(df20) <- c("gbM","TE-like","Unclassified","Unmethylated")
  
  for(f in c('gbM','TE-like','Unclassified','Unmethylated')){
    print(f)
    df21 <- data.frame(table(df11$Classification.x)+table(df11$Classification.y))
    g <- df21[df21$Var1 == f,]$Freq
    if(g < mean(df20[,f])){
      h <- pnorm(g,mean=mean(df20[,f]),sd=sd(df20[,f]),lower.tail=T)*2
      print(h)
    } else {
      h <- pnorm(g,mean=mean(df20[,f]),sd=sd(df20[,f]),lower.tail=F)*2
      print(h)
    }  
  }
  
}




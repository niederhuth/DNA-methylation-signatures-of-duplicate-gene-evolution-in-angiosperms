library(ggplot2)
setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/data/")
setwd("C:/Users/kench/Dropbox/ANALYSIS/GeneDuplication_V2/data") #new laptop

#List species to be analyzed

species = c("Aduranensis","Aipaensis","Alyrata","Athaliana", 
            "Bdistachyon","Boleracea","Bvulgaris","Cclementina","Cpapaya",
            "Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
            "Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
            "Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
            "Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
            "Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
            "Vvinifera","Zmays")

#species = "Athaliana"

KaKs <- Ks <- data.frame()

setcolors <- c("Total"="black","gbM-gbM"="#F94144",
               "gbM-teM"="#F8961E","gbM-unM"="#F9C74F", "unM-unM" = "#90BE60", "unM-teM"="#43AA8B","teM-teM"="#277DA1")


for(a in species){
  #output path
  path1 <- paste("../figures_tables/suppli4/",a,sep="")
  if(!file.exists(path1)){
    dir.create(path1)
  }
  #Duplicate pair methylation similarity
  path2 <- paste("../figures_tables/", a,"/",a, "_Duplicate_pair_met2.csv",sep="")
  #Genes classified by methylation status
  path3 <- paste("../figures_tables/", a,"/",a, "_Kaks_values.csv",sep="")
  
  df1 <- read.csv(path2,header=TRUE) #Dup pairs
  df2 <- read.csv(path3,header=TRUE) #kaks values
  #removing gene pairs where one of the pair had undetermined methylation classification
  df1 <- df1[df1$Similarity != 'Undetermined', ] 
  #merging thew two dataframes to get Dup1-Duplication-Classification-kaksvalues
  df3 <- merge(df1[,c(1,5,7)],df2[,c(2,7,8)],by.x="Duplicate.1",by.y="Feature")
  df3$Classification <- factor(df3$Classification, levels = c("gbM-gbM","gbM-teM","gbM-unM", "unM-unM","unM-teM","teM-teM"))
  df4 <- df3[df3$Duplication != "wgd",]
  df6 <- df3[df3$Duplication == "wgd",]
  
  p <- ggplot(df4) +
    geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
    geom_density(aes(x=Ks,color='Total'),size=1) +
    theme_bw() +
    theme(legend.position = "none") + 
    theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
    theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
    scale_color_manual(values = setcolors)
  ggsave(paste(path1,"/",a,"_SGDs_Ks_fig-test1.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
  
  p <- ggplot(df4) + 
    geom_boxplot(aes(y=Ks,x=Classification), col = c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1"),size=1) + 
    theme_bw()+
    theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
    theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
    theme(axis.title.x = element_blank())+
    scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) + scale_y_continuous(expand=c(0.05,0.05),limits=c(0,5))
  ggsave(paste(path1,"/",a,"_SGDs_Ks_boxplot_test1.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
  
  p <- ggplot(df4) +
    geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
    geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
    theme_bw() +
    theme(legend.position = "none") + 
    theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
    theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
    geom_vline(xintercept=1, linetype="longdash", color="grey55") +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
    scale_color_manual(values = setcolors)
  ggsave(paste(path1,"/",a,"_SGDs_KaKs_fig_test1.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
  
  p <- ggplot(df4) + 
    geom_boxplot(aes(y=Ka.Ks,x=Classification), col = c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1"),size=1) + 
    theme_bw()+
    theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
    theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
    theme(axis.title.x = element_blank())+
    geom_hline(yintercept=1, linetype="longdash", color="grey55") +
    scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) + scale_y_continuous(expand=c(0.05,0.05),limits=c(0,2))
    ggsave(paste(path1,"/",a,"_SGDs_KaKs_boxplot_test1.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
    
    #WGD plots
    
    p <- ggplot(df6) +
      geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
      geom_density(aes(x=Ks,color='Total'),size=1) +
      theme_bw() +
      theme(legend.position = "none") + 
      theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
      theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
      scale_color_manual(values = setcolors)
    ggsave(paste(path1,"/",a,"_WGDs_Ks_fig.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
    
    p <- ggplot(df6) + 
      geom_boxplot(aes(y=Ks,x=Classification), col = c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1"),size=1) + 
      theme_bw()+
      theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
      theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
      theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
      theme(axis.title.x = element_blank())+
      scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) + scale_y_continuous(expand=c(0.05,0.05),limits=c(0,5))
    ggsave(paste(path1,"/",a,"_WGDs_Ks_boxplot.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
    
    p <- ggplot(df6) +
      geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
      geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
      theme_bw() +
      theme(legend.position = "none") + 
      theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
      theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
      geom_vline(xintercept=1, linetype="longdash", color="grey55") +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
      scale_color_manual(values = setcolors)
    ggsave(paste(path1,"/",a,"_WGDs_KaKs_fig.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)
    
    p <- ggplot(df6) + 
      geom_boxplot(aes(y=Ka.Ks,x=Classification), col = c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1"),size=1) + 
      theme_bw()+
      theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
      theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
      theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
      theme(axis.title.x = element_blank())+
      geom_hline(yintercept=1, linetype="longdash", color="grey55") +
      scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) + scale_y_continuous(expand=c(0.05,0.05),limits=c(0,2))
    ggsave(paste(path1,"/",a,"_WGDs_KaKs_boxplot.pdf",sep=""),p,device="pdf", width = 2.5, height = 1.8)  
    
    
}


for(b in c("wgd","proximal","dispersed","tandem","transposed")){
  #Plot the density of genes for each methylation class based Ks
  df5 <- df3[df3$Duplication==b,] 
  p <- ggplot(df5) +
    geom_density(aes(x=Ks,color=Classification),position="dodge",size=1) + 
    geom_density(aes(x=Ks,color='Total'),size=1) +
    theme_bw() +
    #theme(legend.position = "none") + 
    theme(axis.text=element_text(color="black"),
          axis.ticks=element_line(color="black"))+
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous("Ks",expand=c(0,0),limits=c(0,5)) +
    scale_color_manual(values = setcolors)
  ggsave(paste(path1,"/",a,"_",b,"_Ks_fig.pdf",sep=""),p,device="pdf", width = 4, height = 3)
  
  p <- ggplot(df5) + 
    geom_boxplot(aes(y=Ks,x=Classification), col = c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1"),size=1) + 
    theme_bw()+
    theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
    theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
    scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) + scale_y_continuous(expand=c(0,0),limits=c(0,5))
  ggsave(paste(path1,"/",a,"_", b, "_Ks_boxplot.pdf",sep=""),p,device="pdf", width = 4, height = 3)
  
  p <- ggplot(df5) +
    geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
    geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
    theme_bw() +
    #theme(legend.position = "none") + 
    theme(axis.text=element_text(color="black"),
          axis.ticks=element_line(color="black"))+
    geom_vline(xintercept=1, linetype="longdash", color="grey55") +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous("Ka/Ks",expand=c(0,0),limits=c(0,2)) +
    scale_color_manual(values = setcolors)
  ggsave(paste(path1,"/",a,"_",b,"_KaKs_fig.pdf",sep=""),p,device="pdf", width = 4, height = 3)
  
  p <- ggplot(df5) + 
    geom_boxplot(aes(y=Ka.Ks,x=Classification), col = c("#F94144","#F8961E","#F9C74F","#90BE60","#43AA8B","#277DA1"),size=1) + 
    theme_bw()+
    theme(axis.text=element_text(size = 16),axis.ticks=element_line(color="black")) +
    theme(axis.title.x = element_text(size =16),axis.title.y = element_text(size = 16)) +
    theme(axis.text.x = element_text(color = "black", size = 12, angle = 30, hjust = 1))+
    geom_hline(yintercept=1, linetype="longdash", color="grey55") +
    scale_x_discrete(labels=c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")) + scale_y_continuous(expand=c(0,0),limits=c(0,2))
  ggsave(paste(path1,"/",a,"_",b,"_KaKs_boxplot.pdf",sep=""),p,device="pdf", width = 4, height = 3)
    
    
    for(c in c("gbM-gbM","gbM-teM","gbM-unM","unM-unM","unM-teM","teM-teM")){
      if(nrow(df3[df3$Classification == c,]) != 0){
        tmp <- ks.test(df3[df3$Classification == c,]$Ks,
                       sample(df3$Ks,size=nrow(df3[df3$Classification == c,]),replace=FALSE))
        Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
                                  D.statistic=tmp$stat,p.value=tmp$p.value))
        tmp <- ks.test(df3[df3$Classification == c,]$Ka.Ks,
                       sample(df3$Ka.Ks,size=nrow(df3[df3$Classification == c,]),replace=FALSE))
        KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
                                      D.statistic=tmp$stat,p.value=tmp$p.value))
      } else {
        Ks <- rbind(Ks,data.frame(Species=a,Duplication=b,Methylation=c,
                                  D.statistic=NA,p.value=NA))
        KaKs <- rbind(KaKs,data.frame(Species=a,Duplication=b,Methylation=c,
                                      D.statistic=NA,p.value=NA))
      }
    }
  }
  


Ks$p.adjust <- p.adjust(Ks$p.value,method="BH")
KaKs$p.adjust <- p.adjust(KaKs$p.value,method="BH")
#Output K-S results
write.csv(Ks,"../figures_tables/new_Ks-distribution-test.csv",quote=FALSE)
write.csv(KaKs,"../figures_tables/new_KaKs-distribution-test.csv",quote=FALSE)





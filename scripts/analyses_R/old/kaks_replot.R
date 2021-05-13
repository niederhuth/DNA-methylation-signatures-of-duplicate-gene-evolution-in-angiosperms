species = c("A. duranensis","A. ipaensis","A. lyrata","A. thaliana","A. trichopoda",
            "B. distachyon","B. oleracea","B. rapa","B. vulgaris","C. clementina","C. papaya",
            "C. lanatus","C. melo","C. rubella","C. sativus","E. grandis","E. guineensis",
            "E. salsugineum","F. vesca","F. x ananassa","G. max","G. raimondii","L. japonicus",
            "M. acuminata","M. domestica","M. esculenta","M. guttatus","M. truncatula","O. sativa",
            "P. hallii","P. persica","P. trichocarpa","P. virgatum","P. vulgaris","P. x bretschneideri",
            "S. bicolor","S. italica","S. lycopersicum","S. tuberosum","S. viridis","T. cacao",
            "V. vinifera","Z. mays")
df3 <- data.frame()
#for(a in species){
#  df1 <- read.csv(paste(a,"_KaKs_values.csv",sep=""),header=T)
#  df2 <- df1[c(2,3,6,7)]
#  df2$species <- a
#  df2 <- df2[c(5,1,2,3,4)]
#  df3 <- rbind(df3,df2)
#}

for(a in species){
  name <- gsub(" ","",gsub("\\. ","",a))
  df1 <- read.csv(paste(name,"_KaKs_values.csv",sep=""),header=T,row.names=1)
  df1$Duplication <- gsub('wgd','WGD',df1$Duplication)
  df1$Duplication <- gsub('tandem','Tandem',df1$Duplication)
  df1$Duplication <- gsub('transposed','Transposed',df1$Duplication)
  df1$Duplication <- gsub('dispersed','Dispersed',df1$Duplication)
  df1$Duplication <- gsub('proximal','Proximal',df1$Duplication)
  for(b in c('WGD','Transposed','Tandem','Dispersed','Proximal')){
    p <- ggplot(df1[df1$Duplication==b,]) +
      geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
      geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
      geom_vline(xintercept=1, linetype="longdash", color="grey55") +
      theme_bw() +
      theme(axis.text=element_text(color="black"),text=element_text(size=12),legend.position='none',
        axis.ticks=element_line(color="black"),plot.title = element_text(hjust = 0.5))+
      scale_y_continuous(expand=c(0,0)) +
      scale_x_continuous("Ka/Ks",limits=c(0,3)) +
      scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
        "TE-like"="#E69F00","Unmethylated"="#CC79A7","Unclassified"="#999999"),
        breaks=c("Total","gbM","TE-like","Unmethylated","Unclassified")) +
      ggtitle(paste(a,b,sep=" "))
    ggsave(paste("../../dryad/new_KaKs/",name,"_",b,"_KaKs.pdf",sep=""),p,device="pdf",width=2.6,height=2.6)
  }
    p <- ggplot(df1[df1$Duplication!='wgd',]) +
      geom_density(aes(x=Ka.Ks,color=Classification),position="dodge",size=1) + 
      geom_density(aes(x=Ka.Ks,color='Total'),size=1) +
      geom_vline(xintercept=1, linetype="longdash", color="grey55") +
      theme_bw() +
      theme(axis.text=element_text(color="black"),text=element_text(size=12),legend.position='none',
        axis.ticks=element_line(color="black"),plot.title = element_text(hjust = 0.5))+
      scale_y_continuous(expand=c(0,0)) +
      scale_x_continuous("Ka/Ks",limits=c(0,3)) +
      scale_color_manual(values=c("Total"="black","gbM"="#56B4E9",
        "TE-like"="#E69F00","Unmethylated"="#CC79A7","Unclassified"="#999999"),
        breaks=c("Total","gbM","TE-like","Unmethylated","Unclassified")) +
      ggtitle(paste(a,"SGD",sep=" "))
  ggsave(paste("../../dryad/new_KaKs/",name,"_SGD_KaKs.pdf",sep=""),p,device="pdf",width=2.6,height=2.6)
}
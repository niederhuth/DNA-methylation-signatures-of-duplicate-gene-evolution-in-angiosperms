species = c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
            "Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
            "Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
            "Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
            "Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
            "Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
            "Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
            "Vvinifera","Zmays")

FE <- data.frame()
for(a in species){
  df1 <- read.csv(paste("../../figures_tables/SGDs/",a,"_KaKs_values.csv",sep=""),header=T)
  df2 <- data.frame(table(df1[df1$Ka.Ks > 1 & df1$P.Value < 0.05,]$Duplication),table(df1$Duplication))[c(1,2,4)]
  colnames(df2) <- c("Duplication","Positive.Selection","Total")
  df2 <- rbind(df2,data.frame(Duplication=c("Total"),Positive.Selection=sum(df2$Positive.Selection),
                              Total=sum(df2$Total)))
  df2$Percent <- df2$Positive.Selection/df2$Total
  df3 <- data.frame(table(df1[df1$Ka.Ks > 1.1,]$Classification),table(df1$Classification))[c(1,2,4)]
  colnames(df3) <- c("Classification","Positive.Selection","Total")
  df3 <- rbind(df3,data.frame(Classification=c("Total"),Positive.Selection=sum(df3$Positive.Selection),
                              Total=sum(df3$Total)))
  df3$Percent <- df3$Positive.Selection/df3$Total
  df3 <- df3[c(1,2,4,5),]
  df3$order <- c(2,3,4,1)
  p <- ggplot(df3) + 
    geom_bar(aes(y=Percent,x=reorder(Classification,order),
                 fill=Classification),stat="identity") + 
    scale_y_continuous("Percent genes Ka/Ks > 1.1",labels=percent,expand=c(0,0)) + 
    theme_bw() + 
    theme(axis.text=element_text(color="black"),
          axis.ticks=element_line(color="black"),
          legend.position="None") + xlab("") + 
    scale_fill_manual(values=c("Total"="black","gbM"="#56B4E9",
                               "TE-like"="#E69F00","Unmethylated"="#CC79A7"))
  ggsave(paste("../../figures_tables/SGDs/",a,"_positive_selection.pdf",sep=""),p,device="pdf")
  for(i in 1:3){
    FE <- rbind(FE,data.frame(species=c(a),classification=df3[i,1],
                              estimate=fisher.test(matrix(c(df3[i,2],
                                                            df3[i,3]-df3[i,2],
                                                            df3[4,2]-df3[i,2],
                                                            df3[4,3]-df3[4,2]-df5[i,3]),
                                                          nrow=2,ncol=2),alternative="two.sided")$estimate,
                              p.value=fisher.test(matrix(c(df3[i,2],
                                                           df3[i,3]-df3[i,2],
                                                           df3[4,2]-df3[i,2],
                                                           df3[4,3]-df3[4,2]-df3[i,3]),
                                                         nrow=2,ncol=2),alternative="two.sided")$p.value))
  }
}
FE$p.adjust <- p.adjust(FE$p.value,method="BH")

write.csv(FE,"../../figures_tables/SGDs/positive_selection.csv",quote=FALSE,row.names=FALSE)

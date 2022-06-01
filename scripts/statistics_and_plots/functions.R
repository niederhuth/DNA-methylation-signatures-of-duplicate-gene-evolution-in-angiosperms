#Function for making dot plots of Enriched GO terms
GOdotplot <- function(x,fdr=0.05){
  require(ggplot2)
  x=head(x[x$fdr < fdr,],10)
  ggplot(x[x$fdr < fdr,],aes(x=Significant/Annotated,
    y=reorder(Term,Significant/Annotated))) + 
    geom_point(aes(color=fdr,size=Significant)) + 
    theme_bw() +
    scale_color_continuous(low="red",high="blue") +
    xlab("Gene Ratio (Gene List/Annotated Genes)") + 
    ylab("") +
    labs(size="Gene Count",color="FDR") +
    ggtitle("Top 10 Enriched GO Terms")
}

#Function for running topGO on a list of genes
topGO <- function(genelist,goTerms,nodeSize,fdr=0.05,filename,path="",returnData=FALSE){
    require(topGO)
    require(GO.db)
    ifelse(!dir.exists(path),dir.create(path), FALSE)
    BP <- new("topGOdata",description="Biological Process",ontology="BP",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    MF <- new("topGOdata",description="Molecular Function",ontology="MF",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    CC <- new("topGOdata",description="Cellular Compartment",ontology="CC",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    FisherBP <- runTest(BP,algorithm="parentchild",statistic="fisher")
    FisherMF <- runTest(MF,algorithm="parentchild",statistic="fisher")
    FisherCC <- runTest(CC,algorithm="parentchild",statistic="fisher")
    BPgenTable <- GenTable(BP,Fisher=FisherBP,ranksOf="Fisher",
        topNodes=length(score(FisherBP)))
    BPgenTable$Fisher<-gsub("< ","",BPgenTable$Fisher)
    MFgenTable <- GenTable(MF,Fisher=FisherMF,ranksOf="Fisher",
        topNodes=length(score(FisherMF)))
    MFgenTable$Fisher<-gsub("< ","",MFgenTable$Fisher)
    CCgenTable <- GenTable(CC,Fisher=FisherCC,ranksOf="Fisher",
        topNodes=length(score(FisherCC)))
    CCgenTable$Fisher<-gsub("< ","",CCgenTable$Fisher)
    BPgenTable$fdr <- p.adjust(BPgenTable$Fisher,method="BH")
    MFgenTable$fdr <- p.adjust(MFgenTable$Fisher,method="BH")
    CCgenTable$fdr <- p.adjust(CCgenTable$Fisher,method="BH")
    write.csv(BPgenTable[BPgenTable$fdr <= fdr,],paste(path,filename,"_BP.csv",sep=""),
        row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_BP.pdf",sep=""),plot=GOdotplot(BPgenTable,fdr=fdr))
    write.csv(MFgenTable[MFgenTable$fdr <= fdr,],paste(path,filename,"_MF.csv",sep=""),
        row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_MF.pdf",sep=""),plot=GOdotplot(MFgenTable,fdr=fdr))
    write.csv(CCgenTable[CCgenTable$fdr <= fdr,],paste(path,filename,"_CC.csv",sep=""),
        row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_CC.pdf",sep=""),plot=GOdotplot(CCgenTable,fdr=fdr))
    if(returnData){
      return(list(BP=BPgenTable,MF=MFgenTable,CC=CCgenTable))
    }
}






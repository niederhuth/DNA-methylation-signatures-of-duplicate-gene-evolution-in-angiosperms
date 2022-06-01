setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/data")

library(ggplot2)
library(reshape2)
library(scales)

species = "Athaliana"

species = c("Aduranensis","Aipaensis","Alyrata","Athaliana", "Atrichopoda",
            "Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
            "Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
            "Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
            "Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
            "Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
            "Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
            "Vvinifera","Zmays")

for( a in species){
  saving <- paste("../figures_tables/",a,sep="")
  if(!file.exists(saving)){
    dir.create(saving)
  }
 dup_path  <- paste(a,"/dupgen/results-unique/classified_genes.tsv",sep="")
 met_path <- paste(a,"/methylpy/results/", a,"_classified_genes.tsv",sep="")
 met_path2 <- paste(a,"/methylpy/results/",sep="")

df1 <- read.table(met_path,header=TRUE,sep="\t") #methylation data
df2 <- read.table(dup_path,header=TRUE,sep="\t") #Dupgene data

#df1$Feature <- gsub("\\..*","",df1[,1])  # this removes anything after the first decimal in df1 name.
#df1$Feature <- gsub("\\.v3.2*","",df1[,1]) # for individual changes to names
df3 <- merge(df2,df1[,c(1,30)],by.x="Feature",by.y="Feature") #merging the two dataframe (methylation data from column 1 and 23)
write.csv(df3, file =(paste(saving,"/", a, "_MetDup.csv", sep="")))
write.csv(df1, file =(paste(met_path2,"/", a, "_MetClassified_genes.csv", sep="")))

df4 <- as.data.frame(table(df3[2:3])) #function table on column 2 and 3 to count up categories
excel <- as.data.frame(table(df3[2:3]))
excel$Species <- a
#write.csv(excel, file =(paste(saving,"/", a, "_MetDup_Freq.csv", sep="")))

df4$metperc <- df4$Perc2 <- df4$Perc <- NA  #creating three columns with NA as values

for(i in df4$Duplication){
  x = sum(df4[df4$Duplication==i,]$Freq)
  df4[df4$Duplication==i,]$Perc = df4[df4$Duplication==i,]$Freq/x
} #loop to calculate percentage of each duplication and save in the column df4$Perc 

for(i in df4$Classification){
  x = sum(df4[df4$Classification==i,]$Freq)
  df4[df4$Classification==i,]$Perc2 = df4[df4$Classification==i,]$Freq/x
} #loop to calculate percentage of each methylation classification and save in the column df4$Perc2 

for (i in df4$Classification){
  x = sum(df4[df4$Classification==i,]$Freq)
  y = sum(df4[df4$Classification]$Freq)
  df4[df4$Classification==i,]$metperc = x/y
} 
write.csv(df4, file =(paste(saving,"/", a, "_MetDup_Freq.csv", sep="")))

df4$Classification <- factor(df4$Classification, levels = c("gbM", "teM", "unM", "Unclassified"))
df4$Duplication <- factor(df4$Duplication, levels = c("wgd", "tandem", "proximal", "transposed", "dispersed", "singletons", "unclassified"))

p <- ggplot(df4) + 
  geom_bar(aes(x=Duplication,y=Perc,fill=Classification),position="dodge",stat="identity") +
  theme_bw() +
  scale_fill_manual(values=c("#56B4E9", "#E69F00", "#CC79A7", "#999999"))  +  
  #scale_fill_discrete("Methylation Classification") +
  geom_hline(data = df4, aes(yintercept = metperc, color = Classification), linetype="dashed", size = 0.65, show.legend = FALSE) +
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#CC79A7", "#999999")) +
  labs(x="Gene Duplication Type", y= "Percent Genes") +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0))
ggsave(paste(saving,"/",a,"_DupClass.pdf",sep=""),p, device="pdf", width = 8, height = 6)


#Repeat, this time using methylation classification as x-axis, Perc2 as y-axis, and dupgen classification as fill


q <- ggplot(df4) + 
  geom_bar(aes(x=Classification,y=Perc2,fill=Duplication),position="dodge",stat="identity") +
  theme_bw() +
  scale_fill_discrete("Duplication classification") +
  labs(x="Methylation classification", y="Percent Genes") +
  theme(axis.text.x = element_text(color = "black", size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14)) +
  scale_y_continuous(expand = c(0,0))
ggsave(paste(saving,"/",a,"_MetClass.pdf",sep=""),q,device="pdf", width = 8, height = 6)

}

####combining freq datafiles####

setwd("~/Dropbox/ANALYSIS/GeneDuplication_V2/figures_tables/")

species = c("Aduranensis","Aipaensis","Alyrata","Athaliana", "Atrichopoda",
            "Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
            "Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
            "Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
            "Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
            "Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
            "Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
            "Vvinifera","Zmays")

for (a in species){
  path <- paste(a,"/", a,"_MetDup_Freq.csv",sep="")
  path2 <- paste(a,"/", sep="")
  df1 <- read.csv(path, header=TRUE)
  df1 <- cbind(df1, a)
  write.csv(df1, file =paste(path2, a, "_MetDup_Freq2.csv", sep=""))
  } #adding a column with species name to all metDup_freq files 


df2 <- read.csv("./Athaliana/Athaliana_MetDup_Freq2.csv", header=TRUE)

for (a in species){
  path <- paste(a,"/", a,"_MetDup_Freq2.csv",sep="")
  df1 <- read.csv(path, header=TRUE)
  df2 <- rbind(df2, df1)
}

write.csv(df2, file ="All_MetDup_Freq2.csv")


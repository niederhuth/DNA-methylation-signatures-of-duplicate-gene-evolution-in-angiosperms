

#### Creating table for statistical tests ####
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
  df2 <- read.csv(path, header=TRUE)
  
  df2$Col4 <- df2$Col3 <- df2$Col2 <- df2$dupcount <- df2$metcount <- df2$sps <-NA
  total <- sum(df2$Freq)
  df2$sps <- a
  
  for(i in df2$Duplication){
    x = sum(df2[df2$Duplication==i,]$Freq)
    df2[df2$Duplication==i,]$dupcount = x
  } #loop to get each duplication class count
  
  for(i in df2$Classification){
    y = sum(df2[df2$Classification==i,]$Freq)
    df2[df2$Classification==i,]$metcount = y
  } #loop to get each methylation class count
  
  for(i in df2$Duplication){
    b = (df2$dupcount - df2$Freq)
    df2$Col2 = b
  }
  
  for(i in df2$Duplication){
    d = (df2$metcount - df2$Freq)
    df2$Col3 = d
  } 
  
  for(i in df2$Classification){
    c = ((total - df2$dupcount) - (df2$metcount - df2$Freq))
    df2$Col4 = c
  } #creating values for column 4 in the contigency table
  
  write.csv(df2, file =(paste(a,"/", a, "_ForContigencyTable.csv", sep="")))
}

#making one large csv file for contigency table from all species

df3 <- read.csv("All_ContigencyTable.csv") # this could be improved, but for now, I create an empty 
                                           # ouptut file and add all files into this
for (a in species){
  path <- paste(a,"/", a,"_ForContigencyTable.csv",sep="")
  df1 <- read.csv(path, header=TRUE)
  df3 <- rbind(df3, df1)
}
write.csv(df3, file ="All_ContigencyTable.csv")


#### Statistical Test - Fisher-exact test ####

# Individulally testing significant p-value

# test <- matrix(c(1532,	4066,	6334,	17778), nrow = 2,
#                       dimnames =
#                         list(c("gbM", "Not-gbM"),
#                              c("Dispersed", "Not-dispersed")))
# test
# fisher.test(test, alternative = "two.sided")
# fisher.test(Convictions, conf.int = FALSE)
# fisher.test(Convictions, conf.level = 0.95)$conf.int
# fisher.test(Convictions, conf.level = 0.99)$conf.int

#Fisher test for the whole dataset
#from the all_ContigencyTable
#For this a .txt file was created from All_ContigencyTable.csv with just the four columns

alltables <- read.table("ContigencyTable.txt", header = TRUE)
data <- apply(alltables,1, function(x) fisher.test(matrix(x,nr=2), alternative ="two.sided")$p.value)
write.csv(data, file ="MetDup_Fisher_pvalue.csv") 

data <- apply(alltables,1, function(x) fisher.test(matrix(x,nr=2), alternative ="two.sided")$estimate) #getting odds ratio
write.csv(data, file ="MetDup_Fisher_estimate.csv")

data <- apply(alltables,1, function(x) fisher.test(matrix(x,nr=2), alternative ="two.sided")$conf.int) #getting confidence interval ratio
write.csv(data, file ="MetDup_Fisher_confint.csv")


#Multiple test correction using BH

q <- read.csv("MetDup_Fisher_pvalue.csv")
p <- q[,2]
BH <- p.adjust(p, method = "BH", n = length(p))

tmp <- read.csv("MetDup_Fisher_pvalue.csv")
tmp <- cbind(tmp, BH)
write.csv(tmp, "MetDup_Fisher_pvalue_BH.csv")


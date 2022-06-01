library(treeio)
library(ggtree)
library(ape)

species <- c("Aduranensis","Aipaensis","Alyrata","Athaliana","Atrichopoda",
	"Bdistachyon","Boleracea","Brapa","Bvulgaris","Cclementina","Cpapaya",
	"Clanatus","Cmelo","Crubella","Csativus","Egrandis","Eguineensis",
	"Esalsugineum","Fvesca","Fxananassa","Gmax","Graimondii","Ljaponicus",
	"Macuminata","Mdomestica","Mesculenta","Mguttatus","Mtruncatula","Osativa",
	"Phallii","Ppersica","Ptrichocarpa","Pvirgatum","Pvulgaris","Pxbretschneideri",
	"Sbicolor","Sitalica","Slycopersicum","Stuberosum","Sviridis","Tcacao",
	"Vvinifera","Zmays")

tre <- read.tree("orthofinder/orthofinder/Results_Oct03/Species_Tree/SpeciesTree_rooted_node_labels.txt")
tre2 <- as_tibble(tre)

df1 <- read.table("orthofinder/orthofinder/Results_Oct03/Gene_Duplication_Events/Duplications.tsv",header=TRUE,sep="\t")

speciesLevels <- data.frame()
for(i in species){
	a <- b <- tre2[tre2$label==i,]$node
	while(length(a) > 0){
		a <- parent(tre2,a)$node
		b <- gsub("$",paste(",",a,sep=""),b)
	}
	speciesLevels <- rbind(speciesLevels,data.frame(Species=i,Node_Levels=b))
}




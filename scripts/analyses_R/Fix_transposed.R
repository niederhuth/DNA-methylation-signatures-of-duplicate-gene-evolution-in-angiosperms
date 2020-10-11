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

path1 <- paste("orthofinder/orthofinder/",dir("orthofinder/orthofinder/"),sep="")
path2 <- paste(path1,"/Resolved_Gene_Trees/",sep="")
speciesTree <- read.tree(paste(path1,"/Species_Tree/SpeciesTree_rooted_node_labels.txt",sep=""))

geneCount <- data.frame()
for(a in species){
	df1 <- read.table("Athaliana/dupgen/results-unique/Athaliana.transposed.pairs-unique",
		sep="\t",header=TRUE)
	df2 <- read.table("Athaliana/ref/mcscanx/Athaliana_orthogroups.tsv",sep="\t",head=FALSE)

	for(x in c(1:nrow(df1))){
		b <- as.vector(df1[x,]$Transposed)
		p <- df1[x,]$Parent
		og <- df2[df2$V1==b,2]
		tre <- read.tree(paste(path2,og,"_tree.txt",sep=""))
		tre2 <- as_tibble(tre)
		ancestor <- getMRCA(tre,c(paste(a,"_",b,sep=""),paste(a,"_",p,sep="")))
		for(i in c(1:length(child(tre,ancestor)))){
			if(!(tre2[tre2$node==child(tre,ancestor)[i],]$label %in% 
				c(paste(a,"_",b,sep=""),paste(a,"_",p,sep="")))){
				x <- as_tibble(tree_subset(tre,tre2[tre2$node==child(tre,ancestor)[i],]$label,
				levels_back=0))
				if(paste(a,"_",b,sep="") %in% tre2$label){
					c1 <- nrow(x[grep(a,x$label),])
					c2 <- nrow(x[grep("Atrichopoda",x$label),])
					c3 <- nrow(x[grep("Nnucifera",x$label),])
					c4 <- nrow(x[grep("Acoerulea",x$label),])
					c5 <- nrow(x[grep("Spolyrhiza",x$label),])
					tmp <- tree_subset(tre,paste(a,"_",b,sep=""),levels_back=1)
					closest <- tmp$tip.label[grep(paste(a,"_",b,sep=""),
						tmp$tip.label,invert=TRUE)]
				}
			} else if(tre2[tre2$node==child(tre,ancestor)[i],]$label == paste(a,"_",b,sep="")){
				c1 <- c2 <- c3 <- c4 <- c5 <- 0
				closest <- NA
			}
		}
		geneCount <- rbind(geneCount,data.frame(Transposed=b,Parent=p,speciesCount=c1,
			Atrichopoda=c2,Nnucifera=c3,Acoerulea=c4,Spolyrhiza=c5,Closest=closest))
	}
}

#test <- drop.tip(tre,tre$tip.label[!(tre$tip.label %in% x6$label)])


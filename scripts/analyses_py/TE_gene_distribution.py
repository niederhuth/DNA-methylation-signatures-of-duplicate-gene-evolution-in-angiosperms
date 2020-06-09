import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

from functions import chr_filter

#Set variables 
genome_file='ref/'+sys.argv[1]+'.fa.fai'
TEs='ref/annotations/'+sys.argv[1]+'-TEanno.gff'
genes='ref/mcscanx/'+sys.argv[1]+'.gff'
mC_class='methylpy/results/'+sys.argv[1]+'_classified_genes.tsv'
filter_chr=['ChrL','ChrC','ChrM']
window_size=100000
stepsize=50000
output='methylpy/results'+sys.argv[1]+'_TE_gene_distribution.tsv'

#Set chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Create bed file of sliding windows
w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,w=window_size,s=stepsize,i='srcwinnum').filter(chr_filter,chrs)

#Read in genes & coordinates from mcscanx gff
a = pd.read_csv(genes,header=None,sep="\t").replace(sys.argv[1]+'_','',regex=True,inplace=False)
#Read in genes classified by methylation
b = pd.read_csv(mC_class,header=0,sep="\t")
#Merge the two and create a new bedfile
g_bed = pbt.BedTool.from_dataframe(pd.merge(a,b,left_on=1,right_on='Feature')[[0,2,3,1,'Classification']])
#Get the intersections between windows and genes
gene_intersections = pbt.bedtool.BedTool.intersect(w_bed,g_bed,wao=True)
#Get the intersections between windows and TEs
TE_intersections = pbt.bedtool.BedTool.intersect(w_bed,TEs,wao=True)
#Read in gene intersections as a table, keeping only window name & classification
c = pd.read_csv(gene_intersections.fn,header=None,usecols=[3,8],sep="\t")
#Read in TE intersections as a table, keeping only window name and number of nucleotides
d = pd.read_csv(TE_intersections.fn,header=None,usecols=[3,13],sep="\t")
#Create a list of windows
w_list = pd.read_csv(w_bed.fn,header=None,usecols=[3],sep="\t").drop_duplicates().values.tolist()
#List of column names for new table
columns = ['Window','Genes','TEs','TE-nucleotidess','gbM','TE-like','Unmethylated','Unclassified']
#Create empty dataframe with column names
e = pd.DataFrame(columns=columns)
#Iterate over windows list
for i in w_list:
	#For each window, get the number of genes
	Genes = len(c[(c[3]==i[0])])
	#Count number of TEs
	TEs = len(d[d[3]==i[0]])
	#Sum up number of nucleotides that are TE
	TEnts = d[d[3]==i[0]][13].sum()
	#Count number of gbM genes
	gbM = len(c[(c[3]==i[0]) & (c[8]=='gbM')])
	#Count number of TE-like genes
	TElike = len(c[(c[3]==i[0]) & (c[8]=='TE-like')])
	#Count number of Unmethylated genes
	Unmethylated = len(c[(f[3]==i[0]) & (c[8]=='Unmethylated')])
	#Count number of Unclassified genes
	Unclassified = len(c[(f[3]==i[0]) & (c[8]=='Unclassified')])
	#Append to dataframe
	e = e.append(pd.DataFrame([[i[0],Genes,TEs,TEnts,gbM,TElike,Unmethylated,Unclassified]],columns=columns),ignore_index=True)
#Save the results
e.to_csv(output,sep='\t',index=False)





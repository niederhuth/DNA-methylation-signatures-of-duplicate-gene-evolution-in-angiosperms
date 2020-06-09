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
genes="ref/mcscanx/Athaliana.gff"
mC_class='methylpy/results/'+sys.argv[1]+'_classified_genes.tsv'
filter_chr=['ChrL','ChrC','ChrM']
window_size=100000
stepsize=50000


#Set chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Create bed file of sliding windows
w_bed = pbt.bedtool.BedTool.window_maker(pbt.BedTool(genome_file),g=genome_file,w=window_size,s=stepsize,i='srcwinnum').filter(chr_filter,chrs)

#

b = pd.read_csv(genes,header=None,sep="\t").replace(species+'_','',regex=True,inplace=False)
c = pd.read_csv(mC_class,header=0,sep="\t")
d = pbt.BedTool.from_dataframe(pd.merge(b,c,left_on=1,right_on='Feature')[[0,2,3,1,'Classification']])

e = pbt.bedtool.BedTool.intersect(w_bed,d,wao=True)
j = pbt.bedtool.BedTool.intersect(w_bed,TEs,wao=True)
f = pd.read_csv(e.fn,header=None,usecols=[3,8],sep="\t")
k = pd.read_csv(j.fn,header=None,usecols=[3,13],sep="\t")

a = pd.read_csv(w_bed.fn,header=None,usecols=[3],sep="\t").drop_duplicates().values.tolist()

l = ['Window','Genes','TE-nts','gbM','TE-like','Unmethylated','Unclassified']

g = pd.DataFrame(columns=l)
for i in a:
	Genes = len(f[(f[3]==i[0])])
	TEnts = k[k[3]==i[0]][13].sum()
	gbM = len(f[(f[3]==i[0]) & (f[8]=='gbM')])
	TElike = len(f[(f[3]==i[0]) & (f[8]=='TE-like')])
	Unmethylated = len(f[(f[3]==i[0]) & (f[8]=='Unmethylated')])
	Unclassified = len(f[(f[3]==i[0]) & (f[8]=='Unclassified')])
	h = [i[0],Genes,TEnts,gbM,TElike,Unmethylated,Unclassified]
	g = g.append(pd.DataFrame([h],columns=l),ignore_index=True)

g.to_csv('test.tsv', sep='\t', index=False)





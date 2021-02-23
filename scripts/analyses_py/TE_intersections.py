import os
import sys
import pandas as pd
import pybedtools as pbt
import gffutils

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#Set Variables
updown_stream=[100,250,500,750,1000,1500,2000,3000,4000,5000,7500,10000]
annotations='ref/annotations/'+sys.argv[1]+'.gff'
TEs='ref/annotations/'+sys.argv[1]+'-TEanno.gff'
genome_file='ref/'+sys.argv[1]+'.fa.fai'
feature='gene'
filter_chr=['ChrL','ChrC','ChrM']
output='methylpy/results/'+sys.argv[1]+'_TE_intersections.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Create temporary gff database
db = gffutils.create_db(annotations, dbfn='temp.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
pf_bed = ''
for pf in db.features_of_type('gene'):
	pf_bed = pf_bed + '\t'.join([str(pf.seqid), str(pf.start), str(pf.stop), str(pf.id),'.',str(pf.strand)]) + '\n'
	#pf_bed = pf_bed + '\t'.join([str(pf.seqid), str(pf.start), str(pf.stop), str(pf.attributes['Name'][0]),'.',str(pf.strand)]) + '\n'

#
g_bed = pbt.BedTool(pf_bed, from_string=True).filter(functions.chr_filter,chrs).saveas('g_bed.tmp')

TE_table = pd.read_csv(g_bed.fn,header=None,usecols=[3],sep="\t")
TE_table.rename(columns={3:'Gene'},inplace=True)

b='gene_body'
mapping = pbt.bedtool.BedTool.intersect(g_bed,TEs,wa=True,wb=True)
m = pd.read_csv(mapping.fn,header=None,usecols=[3],sep="\t")
tmp = m[3].value_counts().to_dict()
TE_table[b] = TE_table['Gene'].map(tmp)

for a in updown_stream:
	b='up_'+str(a)+'bp'
	s_bed = pbt.bedtool.BedTool.flank(g_bed,g=genome_file,l=a,r=0,s=True).saveas('s_bed.tmp')
	mapping = pbt.bedtool.BedTool.intersect(s_bed,TEs,wa=True,wb=True)
	m = pd.read_csv(mapping.fn,header=None,usecols=[3],sep="\t")
	tmp = m[3].value_counts().to_dict()
	TE_table[b] = TE_table['Gene'].map(tmp)

for a in updown_stream:
	b='down_'+str(a)+'bp'
	s_bed = pbt.bedtool.BedTool.flank(g_bed,g=genome_file,l=0,r=a,s=True).saveas('s_bed.tmp')
	mapping = pbt.bedtool.BedTool.intersect(s_bed,TEs,wa=True,wb=True)
	m = pd.read_csv(mapping.fn,header=None,usecols=[3],sep="\t")
	tmp = m[3].value_counts().to_dict()
	TE_table[b] = TE_table['Gene'].map(tmp)


TE_table.set_index('Gene',inplace=True)
TE_table.to_csv(output, sep='\t', index=True)

tmp=['g_bed.tmp','s_bed.tmp','temp.db']
for b in tmp:
	os.remove(b)










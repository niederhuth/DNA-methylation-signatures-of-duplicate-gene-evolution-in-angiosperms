import os
import sys
import pandas as pd
import pybedtools as pbt

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
annotations='../ref/annotations/'+sys.argv[1]+'.gff'
annotations2='../ref/annotations/'+sys.argv[1]+'_promoter.gff'
genome_file='../ref/'+sys.argv[1]+'.fa.fai'
filter_chr=['ChrL','ChrC','ChrM']
mc_type=['CG','CHG','CHH']
upstream=200
downstream=100
cutoff=3
site_cutoff_only=True
first_feature='gene'
second_feature='gene'
output='results/promoter_methylation.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Make promoter gff file
a = pbt.bedtool.BedTool.flank(annotations,g=genome_file,l=upstream,r=downstream,
	s=True).saveas(annotations2)

#get gene methylation data
print('Getting gene methylation data')
functions.feature_methylation(allc,annotations2,genome_file,output=output,
	mc_type=mc_type,updown_stream=0,feature=first_feature,
	cutoff=cutoff,chrs=chrs,site_cutoff_only=site_cutoff_only)

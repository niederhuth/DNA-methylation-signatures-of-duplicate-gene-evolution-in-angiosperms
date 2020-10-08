import os
import sys
import pandas as pd
import pybedtools

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
annotations='../ref/annotations/'+sys.argv[1]+'.gff'
annotations='../ref/annotations/'+sys.argv[1]+'_promoter.gff'
genome_file2='../ref/'+sys.argv[1]+'.fa.fai'
filter_chr=['ChrL','ChrC','ChrM']
mc_type=['CG','CHG','CHH']
updown_stream=0
cutoff=3
site_cutoff_only=True
first_feature='gene'
second_feature='gene'
output='results/promoter_methylation.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Make promoter gff file


#get gene methylation data
print('Getting gene methylation data')
functions.feature_methylation(filtered_output,annotations,genome_file,output=output,
	mc_type=mc_type,updown_stream=updown_stream,feature=first_feature,
	cutoff=cutoff,chrs=chrs,site_cutoff_only=site_cutoff_only)

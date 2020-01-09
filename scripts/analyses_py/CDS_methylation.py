import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
annotations='../ref/annotations/'+sys.argv[1]+'.gff'
genome_file='../ref/'+sys.argv[1]+'.fa.fai'
filter_chr=['ChrL','ChrC']
mc_type=['CG','CHG','CHH']
updown_stream=0
cutoff=0
first_feature='gene'
second_feature='CDS'
filtered_output="CDS_filtered_allc.tmp"
output='results/CDS_methylation.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#pre-filter data
print('Filtering allc file')
functions.allc_annotation_filter(allc,annotations,genome_file,output=filtered_output,
	updown_stream=updown_stream,first_feature=first_feature,
	second_feature=second_feature,chrs=chrs)

#get gene methylation data
print('Getting gene methylation data')
functions.feature_methylation(filtered_output,annotations,genome_file,output=output,
	mc_type=mc_type,updown_stream=updown_stream,
	feature=first_feature,cutoff=cutoff,chrs=chrs)

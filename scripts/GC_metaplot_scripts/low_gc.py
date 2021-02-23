import os
import sys
import pandas as pd

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='methylpy/allc_'+sys.argv[1]+'.tsv.gz'
annotations='ref/annotations/'+sys.argv[1]+'_low_GC.gff'
genome_file='ref/'+sys.argv[1]+'.fa.fai'
filter_chr=['ChrL','ChrM','ChrC']
mc_type=['CG','CHG','CHH','CH']
window_number=60
updown_stream=2000
cutoff=0
first_feature='gene'
second_feature='CDS'
filtered_data_output='low_gc_filtered_allc.tmp'
output='methylpy/results/'+sys.argv[1]+'_low_gc_metaplot.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#get gene metaplot data
functions.gene_metaplot(allc,annotations,genome_file,output=output,
	mc_type=mc_type,window_number=window_number,updown_stream=updown_stream,
	cutoff=cutoff,first_feature=first_feature,second_feature=second_feature,
	chrs=chrs,filtered_data_output=filtered_data_output,remove_tmp=True)

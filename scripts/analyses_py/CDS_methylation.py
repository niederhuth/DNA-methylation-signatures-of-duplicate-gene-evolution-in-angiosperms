import os
import sys
import pandas as pd

functionsfile = '../../../scripts/new_functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import new_functions as functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
gff='../ref/annotations/'+sys.argv[1]+'.gff'
genome_file='../ref/'+sys.argv[1]+'.fa.fai'
filter_chr=['ChrL','ChrC','ChrM']
mc_type=['CG','CHG','CHH','CH']
updown_stream=0
cutoff=3
site_cutoff_only=True
primary_feature='gene'
secondary_feature='CDS'
bedfile="allc.bed"
output='results/CDS_methylation.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#Convert allc to bed file
print('Converting allc file to bed file')
functions.allc2bed(allc,bedfile)

#get gene methylation data
print('Getting gene methylation data')
functions.feature_methylation(bedfile,gff,output,mc_type=mc_type,
	primary_feature=primary_feature,secondary_feature=secondary_feature,
	cutoff=cutoff,chrs=chrs,site_cutoff_only=site_cutoff_only)

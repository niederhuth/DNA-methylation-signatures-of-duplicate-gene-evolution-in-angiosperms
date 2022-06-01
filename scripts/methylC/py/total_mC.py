import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
genome_file='../ref/'+sys.argv[1]+'.fa.fai'
filter_chr=['ChrL','ChrC','ChrM']
mc_type=['CG','CHG','CHH']
cutoff=0
output='results/'+sys.argv[1]+'_total_weighted_methylation.tsv'

#get chromosome list
chrs = list(pd.read_csv(genome_file,header=None,usecols=[0],dtype='str',sep="\t")[0])
chrs = list(set(chrs).difference(filter_chr))

#get total weighted mC
functions.total_weighted_mC(allc,output=output,mc_type=mc_type,cutoff=cutoff,chrs=chrs)


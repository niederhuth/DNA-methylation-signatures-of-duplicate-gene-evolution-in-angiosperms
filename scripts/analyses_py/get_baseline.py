import os
import sys
import pandas as pd

functionsfile = '../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
mc_type=['CG','CHG','CHH']
baseline={}
cutoff=10
calc_baseline='True'
min_sites=0
qvalue=0.05
df='tmp'
output='tmp2'

#run binomial test
print('Calculating the baseline')
functions.gene_binom_test(df,output=output,mc_type=mc_type,
	baseline=baseline,cutoff=cutoff,calc_baseline=calc_baseline)


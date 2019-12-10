import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
mc_type=['CG','CHG','CHH']
baseline={}
cutoff=10
calc_baseline='True'
min_sites=0
qvalue=0.05
df='results/all_genes_CDS_methylation.txt'
output='results/binomial_test.txt'
output2='results/classified_genes.txt'

#run binomial test
print('Running the binomial test')
functions.gene_binom_test(df,output=output,mc_type=mc_type,
	baseline=baseline,cutoff=cutoff,calc_baseline=calc_baseline)
print('Classifying genes')
functions.classify_genes(output,output=output2,min_sites=min_sites,qvalue=qvalue)


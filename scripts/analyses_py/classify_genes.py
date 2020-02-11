import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
mc_type=['CG','CHG','CHH']
baseline={'CG':0.31490364762044254,'CHG':0.10840507446858585,'CHH':0.06220668638549682}
calc_baseline=False
min_sites=20
qvalue=0.05
df='results/CDS_methylation.tsv'
output='results/binomial_test.tsv'
output2='results/'+sys.argv[1]+'_classified_genes.tsv'

#run binomial test
print('Running the binomial test')
functions.gene_binom_test(df,output=output,mc_type=mc_type,
	baseline=baseline,min_sites=min_sites,calc_baseline=calc_baseline)
print('Classifying genes')
functions.classify_genes(output,output=output2,min_sites=min_sites,qvalue=qvalue)


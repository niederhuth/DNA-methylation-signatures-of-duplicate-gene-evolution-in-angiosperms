import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
mc_type=['CG','CHG','CHH','CH']
baseline={'CG':0.2829252470136191,'CHG':0.05297463403323828,'CHH':0.010693274279710605,'CH':0.020270110077589717}
calc_baseline=False
use_subCH=True
use_CH=True
min_sites=15
qvalue=0.05
uM_cutoff=1
uM_weighted_mC_cutoff=0.02
df='results/CDS_methylation_CH.tsv'
output='results/binomial_test_CH.tsv'
output2='results/'+sys.argv[1]+'_classified_genes_CH.tsv'

#run binomial test
print('Running the binomial test')
functions.gene_binom_test(df,output=output,mc_type=mc_type,
	baseline=baseline,min_sites=min_sites,calc_baseline=calc_baseline)
print('Classifying genes')
functions.classify_genes(output,output=output2,min_sites=min_sites,qvalue=qvalue,uM_cutoff=uM_cutoff,
	uM_weighted_mC_cutoff=uM_weighted_mC_cutoff,use_subCH=use_subCH,use_CH=use_CH)


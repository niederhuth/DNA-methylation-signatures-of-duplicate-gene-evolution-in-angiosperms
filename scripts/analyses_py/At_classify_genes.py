import os
import sys
import pandas as pd

functionsfile = '../../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
mc_type=['CG','CHG','CHH']
baseline={'CG':0.28348908688859714,'CHG':0.055574290523518824,'CHH':0.011431120569854767}
calc_baseline=False
min_sites=20
qvalue=0.05
uM_cutoff=1
uM_weighted_mC_cutoff=0.02
df='CDS_methylation.tsv'
output='binomial_test.tsv'
output2=sys.argv[1]+'_classified_genes.tsv'

#run binomial test
print('Running the binomial test')
functions.gene_binom_test(df,output=output,mc_type=mc_type,
	baseline=baseline,min_sites=min_sites,calc_baseline=calc_baseline)
print('Classifying genes')
functions.classify_genes(output,output=output2,min_sites=min_sites,qvalue=qvalue,
	uM_cutoff=uM_cutoff,uM_weighted_mC_cutoff=uM_weighted_mC_cutoff)


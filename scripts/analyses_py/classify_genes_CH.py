import os
import sys
import pandas as pd

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
mc_type=['CG','CHG','CHH','CH']
baseline={'CG':0.29599914992309934,'CHG':0.046904099965876946,'CHH':0.010190688580620957,'CH':0.01839807637644858}
calc_baseline=False
use_subCH=True
use_CH=True
min_sites=20
qvalue=0.05
uM_cutoff=1
uM_weighted_mC_cutoff=0.02
df='results/CDS_methylation.tsv'
output='results/binomial_test_CH.tsv'
output2='results/'+sys.argv[1]+'_classified_genes_CH.tsv'

#run binomial test
print('Running the binomial test')
functions.gene_binom_test(df,output=output,mc_type=mc_type,
	baseline=baseline,min_sites=min_sites,calc_baseline=calc_baseline)
print('Classifying genes')
functions.classify_genes(output,output=output2,min_sites=min_sites,qvalue=qvalue,uM_cutoff=uM_cutoff,
	uM_weighted_mC_cutoff=uM_weighted_mC_cutoff,use_subCH=use_subCH,use_CH=use_CH)


import os
import sys
import pandas as pd

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
fasta="ref/mcscanx/"+sys.argv[1]+'-cds.fa'
output='ref/mcscanx/'+sys.argv[1]+'-gc123.tsv'

#GC content
print('Getting GC content for genes')
functions.gc123(fasta,output=output)


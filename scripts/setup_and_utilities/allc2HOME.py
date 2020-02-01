import os
import sys

functionsfile = '../../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

#define variables
allc='allc_'+sys.argv[1]+'.tsv.gz'
output=sys.argv[1]+'_allc2HOME.tsv'

#get total weighted mC
functions.allc2HOME(allc,output=output)


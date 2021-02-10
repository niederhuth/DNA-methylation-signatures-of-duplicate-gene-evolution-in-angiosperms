import os
import sys
import pandas as pd

functionsfile = '../../scripts/functions.py'
sys.path.append(os.path.dirname(os.path.expanduser(functionsfile)))

import functions

functions.add_CH('methylpy/results/CDS_methylation.tsv','methylpy/results/CDS_methylation_CH.tsv')

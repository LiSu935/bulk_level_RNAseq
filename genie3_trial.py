import os
import sys

sys.path.append(os.path.abspath("/storage/htc/joshilab/Su_Li/tools_related/genie3"))

from GENIE3 import *

# load data and gene names:


regulators = ['CD19', 'CDH17', 'RAD51', 'OSR2', 'TBX3']
VIM2 = GENIE3(data, gene_names = gene_names, regulators = regulators)
get_link_list(VIM2,gene_names=gene_names,regulators = regulators)

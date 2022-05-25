import os
import sys

sys.path.append(os.path.abspath("/storage/htc/joshilab/Su_Li/tools_related/genie3"))

from GENIE3 import *


tf_file = "/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/Gmax_TF/Gma_TF_list.txt"

tf_list = []

with open(tf_file, 'r') as f:
  for line in f:
    if line.startswith("TF_ID"):
      continue
    else:
      tf = line.rstrip().split('\t')[1]
      tf_list.append(tf)
      
tf_list = list(set(tf_list))  

# load data and gene names:
input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/CortexPhloemXylem_Counts/'
for file in os.listdir(input_dir):
  if file.endswith('_edited.txt'):

    data = loadtxt(input_dir+file, skiprows = 1)
    f = open(input_dir+file)
    gene_names = f.readline()
    f.close()
    gene_names = gene_names.rstrip('\n').split('\t')

    regulators = tf_list
    VIM2 = GENIE3(data, gene_names = gene_names, regulators = regulators)
    get_link_list(VIM2,gene_names=gene_names,regulators = regulators, file_name = input_dir + file.split('.')[0]+'_ranking.txt')

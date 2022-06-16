# This is for only keeping target genes that are DEGs in the original networks, Dr. Xu suggested that TF might be very low or change very little in different tissues,
# so filtering the target genes are more reasonable. 

import pandas as pd
import os
import glob

fi_name = 'Xylem_72h_Bj_edited_GENIE3linkList_0.005.txt'
input_dir = 'C:\\Users\\lsxgf\\OneDrive - University of Missouri\\Gary\\Jaehyo\\CortexPhloemXylem_EdgeR_DEGs\\'

# for lewis:
input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/results/'
deg_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/degs/'
output_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/tf_of_interest/'


def parsing_results(input_dir, fi_name):

    prefix = '_'.join(fi_name.split('_')[:3])
    deg_prefix = '_'.join(fi_name.split('_')[:2])
    DEG_file = glob.glob(deg_dir+deg_prefix+"*.txt")[0]
    deg_list = list(pd.read_table(deg_dir+DEG_file, header=0, sep='\t')['gene'])

    with open(output_dir+prefix+'_links_containing_degs.txt', 'w') as f:
        with open(input_dir+fi_name, 'r') as f1:
            for line in f:
                if line.startswith('""'):
                    continue
                else:
                    line_list = line.rstrip().split('\t')[1:]
                    if line_list[1][1:-1] in deg_list:
                        new_line = '\t'.join([line_list[0][1:-1], line_list[1][1:-1], line_list[2]])
                        f.write(new_line+'\n')
        f1.close()
    f.close()
    
    
for fi in os.listdir(input_dir):
    if fi.endswith("linkList_0.005.txt"):
        print(fi)
        parsing_results(input_dir, fi)
        print("=========================")





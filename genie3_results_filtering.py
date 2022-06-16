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

    dic_tf = {}
    prefix = '_'.join(fi_name.split('_')[:3])
    deg_prefix = '_'.join(fi_name.split('_')[:2])

    DEG_file = glob.glob(deg_dir+deg_prefix+"*.txt")[0]

    deg_list = list(pd.read_table(DEG_file, header=0, sep='\t')['gene'])

    with open(input_dir+fi_name, 'r') as f:
        for line in f:
            if line.startswith('""'):
                continue
            else:
                line_list = line.rstrip().split('\t')[1:]
                if line_list[0] not in dic_tf.keys():
                    dic_tf[line_list[0][1:-1]] = [line_list[1][1:-1]]
                else:
                    dic_tf[line_list[0][1:-1]].append(line_list[1][1:-1])

    with open(output_dir+prefix+'_TF_of_interest.txt', 'w') as f:
        for key in dic_tf.keys():
            if key in deg_list:
                i = 0
                for item in dic_tf[key]:
                    if item in deg_list:
                        i += 1
                if i > 0:
                    f.write(key+'\n')
    f.close()


for fi in os.listdir(input_dir):
    if fi.endswith("linkList_0.005.txt"):
        print(fi)
        parsing_results(input_dir, fi)
        print("=========================")






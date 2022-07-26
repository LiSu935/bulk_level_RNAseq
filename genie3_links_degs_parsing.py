"""
this is for filtering genie3 links that containing degs and 
summarizing the pairs of multi-one and one-multi
"""

import pandas as pd
import glob
import os

list_tf = list(pd.read_table('/Users/lisu/Documents/DrDongXu/Gary_soybean/Cortex_Phloem_Xylem_Bj_TF/TF_list_planttfdb/'
                             'Gma_TF_list.txt', sep='\t', header=0).iloc[:,1:3]. \
               drop_duplicates(ignore_index=True)['Gene_ID'])

dir_links = '/Users/lisu/OneDrive - University of Missouri/Gary/Jaehyo/genie3_links/containing_degs/'
dir_degs = '/Users/lisu/OneDrive - University of Missouri/Gary/Jaehyo/CortexPhloemXylem_EdgeR_DEGs/'

'Cortex_72h_Bj_links_containing_degs.txt'


def parsing(file):
    prefix = '_'.join(file.split('_')[0:2])
    file_path_deg = glob.glob(dir_degs+prefix+"*_exp_value_stat_sig_anno.txt")[0]
    list_deg = list(pd.read_table(file_path_deg, sep='\t', header=0).iloc[:,:3]. \
               drop_duplicates(ignore_index=True)['gene'])
    with open(dir_links+file, 'r') as f:
        dic_degs_Askeys = {}  # regulator(TF)-target(deg):
        dic_tf_Askeys = {} # regulator(TF)-target(deg):
        for line in f:
            line_list = line.rstrip().split('\t')
            if (line_list[0] in list_tf) & (line_list[1] in list_deg):
                if line_list[1] not in dic_degs_Askeys.keys():
                    dic_degs_Askeys[line_list[1]] = [line_list[0]]
                else:
                    dic_degs_Askeys[line_list[1]].append(line_list[0])
                if line_list[0] not in dic_tf_Askeys.keys():
                    dic_tf_Askeys[line_list[0]] = [line_list[1]]
                else:
                    dic_tf_Askeys[line_list[0]].append(line_list[1])

        with open(dir_links+file.split('.txt')[0]+'_degs_Askyes.txt', 'w') as f1:
            f1.write("degs\tKnown_TF\tNumber_of_TF\n")
            for key1 in dic_degs_Askeys.keys():
                value_str = ', '.join(dic_degs_Askeys[key1])
                f1.write(key1+'\t'+value_str+'\t'+str(len(dic_degs_Askeys[key1]))+'\n')
            f1.close()

        with open(dir_links+file.split('.txt')[0]+'_tf_Askyes.txt', 'w') as f2:
            f2.write("Known_TF\tdegs\tNumber_of_degs\n")
            for key2 in dic_tf_Askeys.keys():
                value_str = ', '.join(dic_tf_Askeys[key2])
                f2.write(key2+'\t'+value_str+'\t'+str(len(dic_tf_Askeys[key2]))+'\n')
            f2.close()

        f.close()


for filename in os.listdir(dir_links):
    if filename.endswith("_Bj_links_containing_degs.txt"):
        parsing(filename)


dic = {}
sample_name_list = []
for filename in os.listdir(dir_links):
    if (filename.endswith("_Bj_links_containing_degs.txt")) & (filename.split('.txt')[0] not in sample_name_list):
        sample_name_list.append(filename.split('_Bj_links_containing_degs.txt')[0])

for samplename in sample_name_list:
    degs_df = pd.read_table(dir_links+samplename+'_Bj_links_containing_degs_degs_Askyes.txt', header=0)
    tf_df = pd.read_table(dir_links+samplename+'_Bj_links_containing_degs_tf_Askyes.txt', header=0)
    # this for getting logFC:
    logfc_df = pd.read_table(glob.glob(dir_degs+samplename+"*_exp_value_stat_sig_anno.txt")[0], header=0) \
        [['gene','logFC']].rename({'gene': 'degs'}, axis=1)
    dic[samplename+'_degsK'] = pd.merge(degs_df, logfc_df, on='degs', how='left')
    dic[samplename+'_tfK'] = tf_df

writer = pd.ExcelWriter(dir_links+'All_tissues_containingDEGs_links_ranking.xlsx', engine='xlsxwriter')

for sheet, frame in dic.items(): # .use .items for python 3.X
    frame.to_excel(writer, sheet_name=sheet)

#critical last step
writer.save()









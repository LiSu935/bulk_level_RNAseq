fi_name = 'Cortex_72h_Bj_edited_GENIE3linkList_0.005.txt'
input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/results/'


dic_tf = {}
with open(input_dir+fi_name, 'r') as f:
    for line in f:
        if line.startswith('""'):
            continue
        else:
            line_list = line.rstrip().split('\t')[1:]
            if line_list[0] not in dic_tf.keys():
                dic_tf[line_list[0]] = [(line_list[1], float(line_list[2]))]
            else:
                dic_tf[line_list[0]].append((line_list[1], float(line_list[2])))

print(len(dic_tf))

len_list = []
for tf in dic_tf.keys():
    len_list.append(len(dic_tf[tf]))

print(sorted(len_list)[0])
print(sorted(len_list)[-1])
print(len_list)

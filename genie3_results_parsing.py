import os 

#fi_name = 'Cortex_72h_Bj_edited_GENIE3linkList_0.005.txt'
input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/results/'

def parsing_results(input_dir, fi_name):
    
    dic_tf = {}
    with open(input_dir+fi_name, 'r') as f:
        for line in f:
            if line.startswith('""'):
                continue
            else:
                line_list = line.rstrip().split('\t')[1:]
                if line_list[0] not in dic_tf.keys():
                    #dic_tf[line_list[0]] = [(line_list[1], float(line_list[2]))]
                    dic_tf[line_list[0]] = [line_list[1][1:-1]]
                else:
                    #dic_tf[line_list[0]].append((line_list[1], float(line_list[2])))
                    dic_tf[line_list[0]].append(line_list[1][1:-1])

    print(len(dic_tf))
    
    sorted_dic = sorted(dic_tf.items(), key=lambda x:len(x))
    for i, j in enumerate(sorted_dic):
        if i <20:
            print(j[0]+'\t'+'\t'.join(j[1]))
        
    print("======== links end =================")
    len_list = []
    for tf in dic_tf.keys():
        len_list.append(len(dic_tf[tf]))

    print(sorted(len_list)[0])
    print(sorted(len_list)[-1])
    print(sorted(len_list))
    
for fi in os.listdir(input_dir):
  if fi.endswith(".txt"):
    print(fi)
    parsing_results(input_dir, fi)
    print("=========================")

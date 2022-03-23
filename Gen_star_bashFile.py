out_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/PhloemXylem7221/star_Out_0323'
input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/PhloemXylem7221/data/StaceyG_29_SOL'
readLength = "100"


bash_tail = 'echo "### Ending at: $(date) ###"'

import os 
sample_list = []
for fi in os.listdir(input_dir):
  if fi.endswith('trimmed.fastq.gz'):
    sample_list.append('_'.join(fi.split("_")[0:3]))
    
sample_list = list(set(sample_list))


for sample in sample_list:
  sample_id = sample.split("_")[0]
  sample_fi_r1 = sample+"_R1_001_trimmed.fastq.gz"
  sample_fi_r2 = sample+"_R2_001_trimmed.fastq.gz"
  commandLine = "STAR --runThreadN 12 --genomeDir /group/xulab/Su_Li/Jaehyo/SoybeanJGI_V4/phytozome/Gmax/Wm82.a4.v1/assembly/genome --readFilesCommand zcat ${input_dir}/"+sample_fi_r1+" ${input_dir}/"+sample_fi_r2+" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "+sample_id+" --sjdbGTFfile /group/xulab/Su_Li/Jaehyo/SoybeanJGI_V4/phytozome/Gmax/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang "+readLength
  bash_header = "#!/bin/bash\n"\
  "#-------------------------------------------------------------------------------\n"\
  "#  SBATCH CONFIG\n"\
  "#-------------------------------------------------------------------------------\n"\
  "## resources\n"\
  "#SBATCH -A xulab\n"\
  "#SBATCH --partition hpc3\n"\
  "#SBATCH --cpus-per-task=1\n"\
  "#SBATCH --mem-per-cpu=48G\n"\
  "#SBATCH --time 2-00:00\n"\
  "## labels and outputs\n"\
  "#SBATCH --job-name=spavg-%j.out\n"\
  "#SBATCH --output="+sample_id+"_star-%j.out  # %j is the unique jobID\n"\
  'echo "### Starting at: $(date) ###"\n'\
  "\n\n"\
  "module load miniconda3\n"\
  "source activate rnaseq\n"\
  "module load star/2.7.9a\n"\
  "input_dir='"+input_dir+"'\n"\
  "cd "+out_dir+"\n\n"

  with open(out_dir+"/run_star_"+sample_id+".sh", "w") as f:
    f.write(bash_header+'\n')
    f.write(commandLine+'\n')
    f.write(bash_tail)
  f.close()
  
    
   

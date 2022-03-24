#module load miniconda3
#source activate rnaseq
#module load samtools/samtools-1.8

#cd /group/xulab/Su_Li/Jaehyo/0421_starOut_V4

#samtools index ./M-3Aligned.sortedByCoord.out.bam

#htseq-count -r pos -f bam -t gene -i Name -q -s no M-3Aligned.sortedByCoord.out.bam /group/xulab/Su_Li/Jaehyo/SoybeanJGI_V4/phytozome/Gmax/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 > /group/xulab/Su_Li/Jaehyo/0421_htseqOut_V4_gene/M-3.count.txt



out_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Cuong_RNAseq/flowering_twoAllele/htseq_Out_0323'
input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Cuong_RNAseq/flowering_twoAllele/star_Out_0323'


bash_tail = 'echo "### Ending at: $(date) ###"'

import os 
sample_list = []
for fi in os.listdir(input_dir):
  if fi.endswith('Aligned.sortedByCoord.out.bam'):
    
    sample_id = fi.split('Aligned.sortedByCoord.out.bam')[0]
    commandLine1 = "cd "+input_dir
    commandLine2 = "samtools index ./"+fi
    commandLine3 = "htseq-count -r pos -f bam -t gene -i Name -q -s no "+fi+" /group/xulab/Su_Li/Jaehyo/SoybeanJGI_V4/phytozome/Gmax/Wm82.a4.v1/annotation/Gmax_508_Wm82.a4.v1.gene_exons.gff3 > "+out_dir+"/"+sample_id+".count.txt"
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
    "#SBATCH --job-name=Cuong-%j.out\n"\
    "#SBATCH --output="+sample_id+"_htseq-%j.out  # %j is the unique jobID\n"\
    'echo "### Starting at: $(date) ###"\n'\
    "\n\n"\
    "module load miniconda3\n"\
    "source activate rnaseq\n"\
    "module load samtools/1.14\n"\
    "\n\n"

    with open(out_dir+"/run_htseq_"+sample_id+".sh", "w") as f:
      f.write(bash_header+'\n')
      f.write(commandLine1+'\n')
      f.write(commandLine2+'\n')
      f.write(commandLine3+'\n')
      f.write(bash_tail+'\n')
    f.close()

# this is a trial for genie3 R version
# Source vignettes: https://bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html

#source("/storage/htc/joshilab/Su_Li/tools_related/genie3/GENIE3.R")
# The GENIE3 group has written it into a Bioconductor package. Use the package. 
library(GENIE3)

expr.matrix = t(read.csv("/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/CortexPhloemXylem_Counts/Cortex_72h_Bj_edited.txt", sep = '\t'))
expr.matrix = expr.matrix[rowSums(expr.matrix[])>0,]

tf_df =  read.csv("/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/Gmax_TF/Gma_TF_list.txt", sep = '\t')
tf_list = unique(tf_df$Gene_ID)

# ============================================================================================ # 
# The following are for building the links between all the genes 
# ============================================================================================ # 

# weightMat = GENIE3(expr.matrix)
#print(dim(weightMat))
#print(weightMat[1:5,1:5])


# ============================================================================================ # 
# Restrict the candidate regulators to a subset of genes
# ============================================================================================ # 
regulators <- intersect(tf_list, rownames(expr.matrix))

set.seed(123)
weightMat <- GENIE3(expr.matrix, regulators=regulators)

# ============================================================================================ # 
# Get the list of the regulatory links
# ============================================================================================ # 
linkList <- getLinkList(weightMat)
print(dim(linkList))
linkList <- getLinkList(weightMat, threshold=0.1)
print(dim(linkList))
class(linkList)

output_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/results/'
write.table(linkList, file=paste(output_dir,"linkList_0.1",".txt",sep=""), row.names = TRUE, col.names = NA, sep="\t")







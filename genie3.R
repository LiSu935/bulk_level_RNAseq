# this is a trial for genie3 R version

source("/storage/htc/joshilab/Su_Li/tools_related/genie3/GENIE3.R")

expr.matrix = read.expr.matrix("/storage/htc/joshilab/Su_Li/tools_related/genie3/data.txt", form="rows.are.samples")
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
regulators <- tf_list
set.seed(123)
weightMat <- GENIE3(exprMatr, regulators=regulators)

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







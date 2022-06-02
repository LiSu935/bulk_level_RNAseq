# this is a trial for genie3 R version
# Source vignettes: https://bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html

#source("/storage/htc/joshilab/Su_Li/tools_related/genie3/GENIE3.R")
# The GENIE3 group has written it into a Bioconductor package. Use the package. 
library(GENIE3)

tf_df =  read.csv("/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/Gmax_TF/Gma_TF_list.txt", sep = '\t')
tf_list = unique(tf_df$Gene_ID)

input_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/CortexPhloemXylem_Counts/'
output_dir = '/storage/htc/joshilab/Su_Li/GaryLab/Jaehyo_RNAseq/genie3_jobs/results/'

# ============================================================================================ # 
# The following are for building the links between all the genes 
# weightMat = GENIE3(expr.matrix)
#print(dim(weightMat))
#print(weightMat[1:5,1:5])
# ============================================================================================ # 



# Write into a function
#' @title Runnning the genie3 program
#' @param expr.matrix: A p(row:gene) x n(column:sample) gene expression matrix 
#' @param fi_name: the filename of the original expr.matrix
# ============================================================================================ # 

run_genie3 = function(expr.matrix, fi_name, link_threshold=0.001){

  # Restrict the candidate regulators to a subset of genes
  regulators <- intersect(tf_list, rownames(expr.matrix))

  set.seed(123)
  weightMat <- GENIE3(expr.matrix, regulators=regulators)
  print(dim(weightMat))

  # Get the list of the regulatory links

  # Getting all the links is too heavy for computation.
  #linkList <- getLinkList(weightMat)
  #print(dim(linkList))


  linkList <- getLinkList(weightMat, threshold=link_threshold)
  print(dim(linkList))

  # write the results
  write.table(linkList, file=paste(output_dir,fi_name,"_GENIE3linkList_", link_threshold,".txt",sep=""), row.names = TRUE, col.names = NA, sep="\t")

}


for (count_fi in list.files(path = input_dir))
{
  if (grepl(".txt", count_fi, fixed=TRUE)){
    fi_name = strsplit(count_fi, ".txt")[[1]][1]
    
    expr.matrix = t(read.csv(paste0(input_dir,count_fi),sep='\t'))
    expr.matrix = expr.matrix[rowSums(expr.matrix[])>0,]
    run_genie3 = function(expr.matrix, fi_name, link_threshold=0.001)
  }

}






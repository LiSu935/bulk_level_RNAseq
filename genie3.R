# this is a trial for genie3 R version

source("/storage/htc/joshilab/Su_Li/tools_related/genie3/GENIE3.R")

expr.matrix = read.expr.matrix("/storage/htc/joshilab/Su_Li/tools_related/genie3/data.txt", form="rows.are.samples")
expr.matrix = expr.matrix[rowSums(expr.matrix[])>0,]
weightMat = GENIE3(expr.matrix)
print(dim(weightMat))
print(weightMat[1:5,1:5])


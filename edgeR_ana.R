# this is to use EdgeR to analyze Cuong's RNA-seq data (flowering, two allels)

#BiocManager::install("edgeR")
library(edgeR)
library(statmod)
library(dplyr)
library("writexl")

out_dir = '/Users/lisu/Documents/DrDongXu/Gary_soybean/Cuong/flowering_twoAllele/edgeR_Out_0405/'
anno =  read.csv('/Users/lisu/Documents/DrDongXu/Course/Jaehyo/SoybeanBM_RNAseq/JGI_V4_update/Gmax_508_Wm82.a4.v1.annotation_info_sorted_uni.txt', sep = '\t')

#data = read.table(file="/Users/lisu/Documents/DrDongXu/Course/Jaehyo/PhloemXylem7221/htseq_Out_0323/merged_counts_final.txt", header = T,sep="\t")
data = read.table(file="/Users/lisu/Documents/DrDongXu/Gary_soybean/Cuong/flowering_twoAllele/htseq_Out_0323/final_merged_count.txt", header = T,sep=",")
rownames(data) = data$gene
data = data[,-1]

sample_list = colnames(data)

group = sapply(sample_list,function(i){
  return(strsplit(i, "[.]")[[1]][[1]])})

targets = data.frame(row.names = sample_list,"group"=group)
targets$group = as.factor(targets$group)


targets = targets %>% mutate(group1 =
             case_when((group == "CN_211") ~ "WT",
                       (group == "CN_213") ~ "AAbb",
                       (group == "CN_205") ~ "aabb",
                       (group == "CN_219") ~ "aaBB",
                       )
)
group = targets$group1
targets$group1 = as.factor(targets$group1)

targets = targets %>% 
  rename(
    group_old = group,
    group = group1
  )
class(targets$group)

design = model.matrix(~0+group)
colnames(design) = levels(targets$group)
y <- DGEList(counts=data,group=group)
keep <- filterByExpr(y, design, min.total.count=10)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
fit = glmQLFit(y,design)

# [1] "aabb" "aaBB" "AAbb" "WT" 
my.contrasts = makeContrasts(
  aabb_vs_WT = aabb-WT,
  aaBB_vs_WT = aaBB-WT,
  AAbb_vs_WT = AAbb-WT,
  aabb_vs_aaBB = aabb-aaBB,
  aabb_vs_AAbb = aabb-AAbb,
  aaBB_vs_AAbb = aaBB-AAbb,
  levels = design)

outPut_all_sig = function(qlf_ob, prefix="Phloem_21d_BvsM", N=10000) {
  
  qlf_Phloem_21d=qlf_ob
  FDR <- p.adjust(qlf_Phloem_21d$table$PValue, method="BH")
  sum(FDR < 0.05)
  
  # Here 7000 is a guess. When use this code later, can set it larger. Need to 
  # do filtering on logFC, FDR
  top <- rownames(topTags(qlf_Phloem_21d,p.value=0.05, sort.by = "logFC",n=N))
  cpm_mat = cpm(y)[top,]
  summary(decideTests(qlf_Phloem_21d, lfc=1,p.value=0.05))
  mat1 = as.matrix(topTags(qlf_Phloem_21d,p.value=0.05, sort.by = "logFC",n=N)[[1]])
  total_mat = cbind(cpm_mat,mat1)
  
  total_mat = cbind(gene = rownames(total_mat), total_mat)
  total_df = as.data.frame(total_mat)
  total_df_anno = total_df %>%
    left_join(y = anno[,-c(1,4)], by = c("gene" = "locusName"))
  
  total_df_anno[,c(2:ncol(total_df))] = sapply(total_df_anno[,c(2:ncol(total_df))], as.numeric)
  #write.table(total_df_anno,file=paste0(out_dir,"Phloem_21d_BvsM_exp_value_stat_anno.txt"),sep = '\t', row.names = F)
  write_xlsx(total_df_anno,paste0(out_dir,prefix,"_exp_value_stat_anno.xlsx"))
  
  total_df_anno_sig = total_df_anno[(total_df_anno$FDR<=0.05)&((total_df_anno$logFC>=1)|(total_df_anno$logFC<=-1)),]
  print(nrow(total_df_anno_sig))
  
  #write.table(total_df_anno_sig,file=paste0(out_dir,"Phloem_21d_BvsM_exp_value_stat_sig_anno.txt"),sep = '\t', row.names = F)
  write_xlsx(total_df_anno_sig,paste0(out_dir,prefix,"_exp_value_stat_sig_anno.xlsx"))
  
}


# aabb_vs_WT = aabb-WT,
#aaBB_vs_WT = aaBB-WT,
#AAbb_vs_WT = AAbb-WT,
#aabb_vs_aaBB = aabb-aaBB,
#aabb_vs_AAbb = aabb-AAbb,
#aaBB_vs_AAbb = aaBB-AAbb,

qlf_aabb_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,"aabb_vs_WT"])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aabb_vs_WT, prefix="aabb_vs_WT")

for (con in c('aabb_vs_WT', "aaBB_vs_WT", "AAbb_vs_WT", "aabb_vs_aaBB", "aabb_vs_AAbb", "aaBB_vs_AAbb")) {
  qlf_ob = glmQLFTest(fit, contrast=my.contrasts[,con])
  prefix_str = paste0("Com_",con)
  outPut_all_sig(qlf_ob, prefix=prefix_str)
}

qlf_aaBB_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,1])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aaBB_vs_WT, prefix="Com_205aabb_vs_211WT")


qlf_aaBB_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,2])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aaBB_vs_WT, prefix="Com_219aaBB_vs_211WT")

qlf_aaBB_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,"AAbb_vs_WT"])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aaBB_vs_WT, prefix="Com_213AAbb_vs_211WT")

qlf_aaBB_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,"aabb_vs_aaBB"])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aaBB_vs_WT, prefix="Com_205aabb_vs_219aaBB")

qlf_aaBB_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,"aabb_vs_AAbb"])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aaBB_vs_WT, prefix="Com_205aabb_vs_213AAbb")

qlf_aaBB_vs_WT = glmQLFTest(fit, contrast=my.contrasts[,"aaBB_vs_AAbb"])
#topTags(qlf_Phloem_21d)
outPut_all_sig(qlf_aaBB_vs_WT, prefix="Com_219aaBB_vs_213AAbb")



# ================================================================================ #
# sessionInfo()
#R version 4.0.4 (2021-02-15)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 10.16

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets 
#[6] methods   base     

#other attached packages:
#[1] writexl_1.4.0  dplyr_1.0.8    statmod_1.4.36
#[4] edgeR_3.32.1   limma_3.46.0  

#loaded via a namespace (and not attached):
 #[1] Rcpp_1.0.8.3        magrittr_2.0.2     
 #[3] tidyselect_1.1.2    lattice_0.20-45    
 #[5] R6_2.5.1            rlang_1.0.2        
 #[7] fansi_1.0.3         tools_4.0.4        
 #[9] grid_4.0.4          xfun_0.30          
#[11] tinytex_0.37        utf8_1.2.2         
#[13] cli_3.2.0           ellipsis_0.3.2     
#[15] tibble_3.1.6        lifecycle_1.0.1    
#[17] crayon_1.5.1        BiocManager_1.30.16
#[19] purrr_0.3.4         vctrs_0.3.8        
#[21] glue_1.6.2          compiler_4.0.4     
#[23] pillar_1.7.0        generics_0.1.2     
#[25] locfit_1.5-9.4      pkgconfig_2.0.3    

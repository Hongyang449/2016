#+ include=FALSE
load('/Users/hyangl/project/ras/results/2016/info/equi_positions.RData')
load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/ras/sca_pairs.RData")
rownames(common_pairs)=NULL 

#+ ras
# ras
common_pairs

#+ include=FALSE
load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/gt/sca_pairs.RData")
rownames(common_pairs)=NULL

#+ gt
# gt
cbind(common_pairs, "a_ras"=gt2ras[common_pairs[,'a']], "b_ras"=gt2ras[common_pairs[,'b']])

#+ include=FALSE
load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/eftu/sca_pairs.RData")
rownames(common_pairs)=NULL 

#+ eftu
# eftu
cbind(common_pairs, "a_ras"=eftu2ras[common_pairs[,'a']], "b_ras"=eftu2ras[common_pairs[,'b']])

#+ close, include=TRUE, eval=FALSE
library(rmarkdown)
render("common_pairs.r", "pdf_document", clean=TRUE)


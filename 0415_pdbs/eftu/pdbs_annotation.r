## name: pdbs_annotation.r
## date: 09/21/2016


load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pdbs.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/eftu_annotate.RData")

tbl_eftu <- eftu_annotate[,c("ligs","colorLig","source")]
colnames(tbl_eftu) <- c("ligandType","ligandColors","source")

save(tbl_eftu, pdbs_eftu_new, eftu_annotate,
     file="pdbs_annotation.RData")




## name: sse_ki.r
## date: 07/16/2016

pdb <- read.pdb("/Users/hyangl/project/ras/results/2016/0601_wt_ki/dccm/pdb/ki.pdb")
sse_kin5 <- dssp(pdb, exefile="mkdssp")
save(sse_kin5, file="sse_ki.RData")



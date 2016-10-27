## name: sse_eftu.r
## date: 07/21/2016

# previous sse can be used for pca
# but for cij we should calculate sse from 1TUI_A (residues 9:405; 397 in total) !!

pdb <- read.pdb("pdbs/split_chain/1TUI_A.pdb")
sse_eftu <- dssp(pdb, exefile="mkdssp")
sse_eftu$helix$start <- sse_eftu$helix$start - 8
sse_eftu$helix$end <- sse_eftu$helix$end - 8
sse_eftu$sheet$start <- sse_eftu$sheet$start - 8
sse_eftu$sheet$end <- sse_eftu$sheet$end - 8

# sse for PL,SI,SII,a3
helix <- list(start=18-8, end=25-8)
helix$length <- helix$end - helix$start + 1
PL <- list(helix=helix)
helix <- list(start=39-8, end=65-8)
helix$length <- helix$end - helix$start + 1
SI <- list(helix=helix)
helix <- list(start=81-8, end=99-8)
helix$length <- helix$end - helix$start + 1
SII <- list(helix=helix)
helix <- list(start=110-8, end=129-8)
helix$length <- helix$end - helix$start + 1
a3 <- list(helix=helix)

# add the concise one - excluding small unimportant sse
sse_eftu_concise <- sse_eftu
inds1 <- -c(2,8,9)
sse_eftu_concise$helix$start <- sse_eftu_concise$helix$start[inds1]
sse_eftu_concise$helix$end <- sse_eftu_concise$helix$end[inds1]
sse_eftu_concise$helix$length <- sse_eftu_concise$helix$length[inds1]
inds2 <- -c(2,3,5:10,12,13,15,17,19,22,23)
sse_eftu_concise$sheet$start <- sse_eftu_concise$sheet$start[inds2]
sse_eftu_concise$sheet$end <- sse_eftu_concise$sheet$end[inds2]
sse_eftu_concise$sheet$length <- sse_eftu_concise$sheet$length[inds2]

# sse_eftu385 for non-gap pdbs positions
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pdbs.RData")
sse_eftu385 <- bounds.sse(pdbs_eftu_new$sse["1TTT_A",!is.gap(pdbs_eftu_new)])
sse_eftu385$sse <- pdbs_eftu_new$sse["1TTT_A",!is.gap(pdbs_eftu_new)]
# add length to sse for draw.box()
sse_eftu385$helix$length <- sse_eftu385$helix$end - sse_eftu385$helix$start + 1
sse_eftu385$sheet$length <- sse_eftu385$sheet$end - sse_eftu385$sheet$start + 1

save(sse_eftu, sse_eftu_concise, sse_eftu385,
     PL, SI, SII, a3,
     file="sse_eftu.RData")



## name: sse_ras.r
## date: 04/18/2016

## name: sse_for_plot.r
## date: 01/21/2015

## Based on the sse of 5P21 (166 residues), I modify it for further plotting

pdb <- read.pdb("pdbs/5P21.pdb")
sse <- dssp(pdb)
sse$sse[92] <- "H"
sse$sse[65:67] <- rep("H",3)

sse$helix$start <- c(16,65,87,127,152)
sse$helix$end <- c(25,74,103,137,164)
sse$helix$length <- c(10,10,17,11,13)
names(sse$helix$start) <- 1:5
names(sse$helix$end) <- 1:5
names(sse$helix$length) <- 1:5
sse$helix$chain <- sse$helix$chain[1:5]
sse$helix$type <- sse$helix$type[1:5]

label_helix <- (sse$helix$start + sse$helix$end)/2
label_sheet <- (sse$sheet$start + sse$sheet$end)/2

names(label_helix) <- paste0("a",names(label_helix))
names(label_sheet) <- paste0("b",names(label_sheet))

label_ab <- c(sse$helix$start, sse$helix$end, sse$sheet$start, sse$sheet$end)

load("/Users/hyangl/project/ras/results/2015/0105_manuscript/pca/sse_for_plot.RData")

sse_ras <- sse

# sse for PL,SI,SII,a3
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
# use which(membership_ras==2)
helix <- list(start=10, end=17)
helix$length <- helix$end - helix$start + 1
PL <- list(helix=helix)
helix <- list(start=25, end=37)
helix$length <- helix$end - helix$start + 1
SI <- list(helix=helix)
helix <- list(start=57, end=75)
helix$length <- helix$end - helix$start + 1
SII <- list(helix=helix)
helix <- list(start=86, end=109)
helix$length <- helix$end - helix$start + 1
a3 <- list(helix=helix)

save(sse_ras,label_helix,label_sheet,label_ab,
     PL,SI,SII,a3,
     file="sse_ras.RData")




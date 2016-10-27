## name: cmap_pdbs.r
## date: 09/09/2016

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs.RData")

source("/Users/hyangl/project/ras/results/2015/functions/cmap.pdbs.R")

pwd <- "/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs/split_chain/"
cutoff_dist <- seq(3,6,by=0.5)

cmap_ras_noh <- lapply(as.list(cutoff_dist), function(x) { 
  cmap.pdbs(pdbs_ras_new, pwd=pwd, dcut=x, scut=2, mask.lower=FALSE)
})
names(cmap_ras_noh) <- cutoff_dist

save(cmap_ras_noh,
     file="cmap_pdbs.RData")



## name: pca_of_pdbs.r
## date: 10/28/2016

## Here I calcuate different distance matrices and compare pca of them.

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pdbs.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pca.RData")
source("/Users/hyangl/project/ras/results/2015/functions/mask.dm.R")

pdbs_eftu_aa <- read.all(pdbs_eftu_new, prefix="/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pdbs/corefit/",
  pdbext=".pdb_flsq.pdb", rm.ligand=TRUE, ncore=NULL)

gaps.res <- gap.inspect(pdbs_eftu_new$ali)
gaps.pos <- gap.inspect(pdbs_eftu_new$xyz)

# 1. Calpha coordinates of core fitted structures
xyz_eftu <- fit.xyz(pdbs_eftu_new$xyz["1TTT_A", ], pdbs_eftu_new, core$c1A.xyz, core$c1A.xyz)
pc_xyz_eftu <- pca(xyz_eftu[, gaps.pos$f.inds])

# 2. Calpha distance matrices (dm)
dmca_eftu <- dm(pdbs_eftu_new)
pc_dmca_eftu <- pca.array(dmca_eftu[gaps.res$f.inds, gaps.res$f.inds, ])

# 3. all-atom distance matrices
dmaa_eftu <- dm(pdbs_eftu_aa$all, grpby=pdbs_eftu_aa$all.grpby, gc.first=TRUE, ncore=10)
pc_dmaa_eftu <- pca.array(dmaa_eftu[gaps.res$f.inds, gaps.res$f.inds, ])

# 4. all-atom dm treated with 10 angstrom cutoff
tdm <- apply(dmaa_eftu, 1:2, function(x) {
  if(all(!is.na(x)) && all(x>10)) {
    rep(0, length(x))
  }
  else {
    x
  }
})
dmaa_cut10_eftu <- apply(tdm, c(3,1), function(x) x); class(dmaa_cut10_eftu) <- 'dmat'
pc_dmaa_cut10_eftu  <- pca.array(dmaa_cut10_eftu[gaps.res$f.inds, gaps.res$f.inds, ])

# 5. all-atom dm treated with a smoothing mask function.
dmaa_masked_eftu <- mask.dm(dmaa_eftu, 4, 8)
pc_dmaa_masked_eftu <- pca.array(dmaa_masked_eftu[gaps.res$f.inds, gaps.res$f.inds, ])

# 6. all-atom cmap
cmaa_eftu <- cmap(pdbs_eftu_aa$all, grpby=pdbs_eftu_aa$all.grpby, dcut=4.5, 
  scut=1, collapse=FALSE, gc.first=TRUE, ncore=10)
pc_cmaa_eftu <- pca.array(cmaa_eftu[gaps.res$f.inds, gaps.res$f.inds,])

# Grouping based on PCA of Cartesian coordinates
hc <- hclust(dist(pc_xyz_eftu$z[, 1:2]))
grps <- cutree(hc, k=2)

layout(matrix(1:6, 2, 3))
plot(pc_xyz_eftu, pc.axes=c(1,2), col=grps, main='Cartesian Coordinates', cex=0.5)
plot(pc_dmca_eftu, pc.axes=c(1,2), col=grps, main='Ca-Ca DM', cex=0.5)
plot(pc_dmaa_eftu, pc.axes=c(1,2), col=grps, main='All-atom DM', cex=0.5)
plot(pc_dmaa_cut10_eftu, pc.axes=c(1,2), col=grps, main='All-atom DM (Cutoff 10 Angstrom)', cex=0.5)
plot(pc_dmaa_masked_eftu, pc.axes=c(1,2), col=grps, main='All-atom DM (Masked 4-8)', cex=0.5)
plot(pc_cmaa_eftu, pc.axes=c(1,2), col=grps, main='Contact Map', cex=0.5)
dev.copy2pdf(file="figures/pca_of_pdbs_eftu.pdf")

dmaa_eftu385 <- dmaa_eftu[gaps.res$f.inds, gaps.res$f.inds,]
# assign names coresponding to 1TTT_A
rownames(dmaa_eftu385) <- pdbs_eftu_new$resno["1TTT_A",gaps.res$f.inds]
colnames(dmaa_eftu385) <- pdbs_eftu_new$resno["1TTT_A",gaps.res$f.inds]

save(pdbs_eftu_aa, pdbs_eftu_new, gaps.res, gaps.pos, grps, dmaa_eftu385,
     xyz_eftu, dmca_eftu, dmaa_eftu, dmaa_cut10_eftu, dmaa_masked_eftu, cmaa_eftu,
     pc_xyz_eftu, pc_dmca_eftu, pc_dmaa_eftu, pc_dmaa_cut10_eftu, pc_dmaa_masked_eftu, pc_cmaa_eftu,
     file="pca_of_pdbs.RData")

save(pdbs_eftu_aa, gaps.res, gaps.pos, grps, dmaa_eftu385, dmaa_eftu, dmaa_masked_eftu,
     file="dmaa.RData")


## name: pca_of_pdbs.r
## date: 10/28/2016

## Here I calcuate different distance matrices and compare pca of them.

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pca.RData")
source("/Users/hyangl/project/ras/results/2015/functions/mask.dm.R")

pdbs_ras_aa <- read.all(pdbs_ras_new, prefix="/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/corefit/",
  pdbext=".pdb_flsq.pdb", rm.ligand=TRUE, ncore=NULL)

gaps.res <- gap.inspect(pdbs_ras_new$ali)
gaps.pos <- gap.inspect(pdbs_ras_new$xyz)

# 1. Calpha coordinates of core fitted structures
xyz_ras <- fit.xyz(pdbs_ras_new$xyz[1, ], pdbs_ras_new, core$c1A.xyz, core$c1A.xyz)
pc_xyz_ras <- pca(xyz_ras[, gaps.pos$f.inds])

# 2. Calpha distance matrices (dm)
dmca_ras <- dm(pdbs_ras_new)
pc_dmca_ras <- pca.array(dmca_ras[gaps.res$f.inds, gaps.res$f.inds, ])

# 3. all-atom distance matrices
dmaa_ras <- dm(pdbs_ras_aa$all, grpby=pdbs_ras_aa$all.grpby, gc.first=TRUE, ncore=10)
pc_dmaa_ras <- pca.array(dmaa_ras[gaps.res$f.inds, gaps.res$f.inds, ])

# 4. all-atom dm treated with 10 angstrom cutoff
tdm <- apply(dmaa_ras, 1:2, function(x) {
  if(all(!is.na(x)) && all(x>10)) {
    rep(0, length(x))
  }
  else {
    x
  }
})
dmaa_cut10_ras <- apply(tdm, c(3,1), function(x) x); class(dmaa_cut10_ras) <- 'dmat'
pc_dmaa_cut10_ras  <- pca.array(dmaa_cut10_ras[gaps.res$f.inds, gaps.res$f.inds, ])

# 5. all-atom dm treated with a smoothing mask function.
dmaa_masked_ras <- mask.dm(dmaa_ras, 4, 8)
pc_dmaa_masked_ras <- pca.array(dmaa_masked_ras[gaps.res$f.inds, gaps.res$f.inds, ])

# 6. all-atom cmap
cmaa_ras <- cmap(pdbs_ras_aa$all, grpby=pdbs_ras_aa$all.grpby, dcut=4.5, 
  scut=1, collapse=FALSE, gc.first=TRUE, ncore=10)
pc_cmaa_ras <- pca.array(cmaa_ras[gaps.res$f.inds, gaps.res$f.inds,])

# Grouping based on PCA of Cartesian coordinates
hc <- hclust(dist(pc_xyz_ras$z[, 1:2]))
grps <- cutree(hc, k=4)

layout(matrix(1:6, 2, 3))
plot(pc_xyz_ras, pc.axes=c(1,2), col=grps, main='Cartesian Coordinates', cex=0.5)
plot(pc_dmca_ras, pc.axes=c(1,2), col=grps, main='Ca-Ca DM', cex=0.5)
plot(pc_dmaa_ras, pc.axes=c(1,2), col=grps, main='All-atom DM', cex=0.5)
plot(pc_dmaa_cut10_ras, pc.axes=c(1,2), col=grps, main='All-atom DM (Cutoff 10 Angstrom)', cex=0.5)
plot(pc_dmaa_masked_ras, pc.axes=c(1,2), col=grps, main='All-atom DM (Masked 4-8)', cex=0.5)
plot(pc_cmaa_ras, pc.axes=c(1,2), col=grps, main='Contact Map', cex=0.5)
dev.copy2pdf(file="figures/pca_of_pdbs_ras.pdf")

dmaa_ras166 <- dmaa_ras[gaps.res$f.inds, gaps.res$f.inds,]

save(pdbs_ras_aa, pdbs_ras_new, gaps.res, gaps.pos, grps, dmaa_ras166,
     xyz_ras, dmca_ras, dmaa_ras, dmaa_cut10_ras, dmaa_masked_ras, cmaa_ras,
     pc_xyz_ras, pc_dmca_ras, pc_dmaa_ras, pc_dmaa_cut10_ras, pc_dmaa_masked_ras, pc_cmaa_ras,
     file="pca_of_pdbs.RData")

save(pdbs_ras_aa, gaps.res, gaps.pos, grps, dmaa_ras166, dmaa_ras, dmaa_masked_ras,
     file="dmaa.RData")


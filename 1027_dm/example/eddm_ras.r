## name: eddm_ras.r
## date: 10/28/2016

#'---
#'title: "Ensemble Difference Distance Matrices analysis of Ras Crystallographic Structures"
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(fig.path='figures/', dev='png', dev.args=list(type="cairo"), dpi=300)
library(bio3d)
library(png)

#+ preamble, include=FALSE, eval=FALSE
library(rmarkdown)
render("eddm_ga3.r", "pdf_document", clean=FALSE)

#+ load pdbs
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs.RData")
source("/Users/hyangl/bio3d/new_funs/eddm_funs.R")

#'
#' ## 1. PCA of crystallographic structures
#'
#' First, try different types of PCA to analyze and group structures. This include conventional
#' Cartesian coordinate PCA with superimposed structures, PCA of C-alpha distance matrix (DMs),
#' PCA of all-atom DM, PCA of all-atom DM with distal pairs filtered by a distance cutoff,
#' PCA of all-atom DM with distal pairs masked with a smooth function, and
#' PCA of contact maps.
#'

gaps.res <- gap.inspect(pdbs_ras_new$ali)
gaps.pos <- gap.inspect(pdbs_ras_new$xyz)

#+ pca_xyz
# PCA of Cartesian coordinates of core fitted structures.
xyz <- fit.xyz(pdbs_ras_new$xyz[1, ], pdbs_ras_new, core$c1A.xyz, core$c1A.xyz)
pc <- pca(xyz[, gaps.pos$f.inds])

#+ pca_dm
# PCA of Calpha-Calpha distance matrices.
dm.ca <- dm(pdbs_ras_new)
pc.ca <- pca.array(dm.ca[gaps.res$f.inds, gaps.res$f.inds, ])

#+ pca_aa_dm
# PCA of residue-wise all-atom distance matrices.
pdbs.aa <- read.all(pdbs_ras_new, prefix="/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/corefit/",
  pdbext=".pdb_flsq.pdb",rm.ligand=TRUE, ncore=NULL)
dm.aa <- dm(pdbs.aa$all, grpby=pdbs.aa$all.grpby, gc.first=TRUE, ncore=10)
pc.aa <- pca.array(dm.aa[gaps.res$f.inds, gaps.res$f.inds, ])

#+ pca_aa_dm_cut (10 angstrom)
# PCA of all-atom DMs treated with a distance cutoff.
tdm <- apply(dm.aa, 1:2, function(x) {
  if(all(!is.na(x)) && all(x>10)) {
    rep(0, length(x))
  }
  else {
    x
  }
})
dm.aa.cut <- apply(tdm, c(3,1), function(x) x)
class(dm.aa.cut) <- 'dmat'
pc.aa.cut <- pca.array(dm.aa.cut[gaps.res$f.inds, gaps.res$f.inds, ])

#+ pca_aa_dm_masked
# PCA of all-atom DMs treated with a smooth masking function.
dm.aa.masked <- mask.dm(dm.aa, 4, 8)
pc.aa.masked <- pca.array(dm.aa.masked[gaps.res$f.inds, gaps.res$f.inds, ])

#+ pca_aa_dm_inverse (for future network building)
dm.aa.inverse <- 1/dm.aa.masked - 1/mask.dm(8, 4, 8)

#+ pca_cmap
# PCA of contact maps.
cm <- cmap(pdbs.aa$all, grpby=pdbs.aa$all.grpby, scut=1, collapse=FALSE, gc.first=TRUE, ncore=10)
pc.cm <- pca.array(cm[gaps.res$f.inds, gaps.res$f.inds,])

#+ plot_cmp_pca, fig.cap='Clustering of structures based on various types of PCA schemes.'
# Grouping based on PCA of Cartesian coordinates
hc <- hclust(dist(pc$z[, 1:2]))
grps <- cutree(hc, k=4)
layout(matrix(1:6, 2, 3))
plot(pc, pc.axes=c(1,2), col=grps, main='Cartesian', cex=0.5)
plot(pc.ca, pc.axes=c(1,2), col=grps, main='Ca-Ca DM', cex=0.5)
plot(pc.aa, pc.axes=c(1,2), col=grps, main='All-atom DM', cex=0.5)
plot(pc.aa.cut, pc.axes=c(1,2), col=grps, main='All-atom DM (Cutoff)', cex=0.5)
plot(pc.aa.masked, pc.axes=c(1,2), col=grps, main='All-atom DM (Mask)', cex=0.5)
plot(pc.cm, pc.axes=c(1,2), col=grps, main='Contact Map', cex=0.5)

#'
#' PC loadings of DM based PCA indicate regions that have the major contribution to the
#' variance of residue distance along each PC. Here, we focus on PCA of "masked" DMs.
#+ plot_aaDM_loading1, fig.cap='PC1 loadings of PCA of (masked) all-atom DM. Red indicate pairs close each other in the "GTP" cluster.'
lm1 <- plot.matrix.loadings(pc.aa.masked, pc=1, plot=FALSE)
gg.mat(lm1) + gg_sse(pdbs$sse[1, gaps.res$f.inds], side=c(1:4)) +
  geom_tile(aes(width=abs(value*5), height=abs(value*5)))

#+ plot_aaDM_loading2, fig.cap='PC2 loadings of PCA of (masked) all-atom DM. Blue indicate pairs deviate each other in the "GDP" cluster.'
lm2 <- plot.matrix.loadings(pc.aa.masked, pc=2, plot=FALSE)
gg.mat(lm2) + gg_sse(pdbs$sse[1, gaps.res$f.inds], side=c(1:4)) +
  geom_tile(aes(width=abs(value*5), height=abs(value*5)))

#'
#' ## 2. Significance of contacts
#'

#'
#' Function **eddm2()** will compute significance of differences in residue contacts
#' between defined groups. Output is a data.frame of contacts data with a row per contacting
#' residue pair and columns containing:
#' (**1**) contact bit where 1 denote a contact between the two residues
#' (columns are called `c1`, `c2`, `c3`, for group 1, 2, 3, respectively),
#' (**2**) mean minimal distances between the pair of residues for each group
#' (columns are denoted `d1`, `d2`, `d3` for groups 1, 2, 3, respectively,
#' (**3**) difference in means between groups
#' (`dm1_2` for the difference between group 1 and 2),
#' (**4**) and p-values
#' (`pv1_2` for the p-value for this contact for group 1 and 2, `pv1_3` for group 1 and 3, etc).
#'

#+ eddm2, cache=TRUE
res2b <- eddm2(pdbs.aa, dm=dm.aa.masked, grps=grps)
head(res2b)

# extract data for grps 1 and 3 only
res3b <- subset.eddm(res2b, grps=c(1,3), alpha=0.005, beta=1.0)
head(res3b)

# add SSE annotations to the table
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
community_ras <- c("a1/b1-b3", "PL", "SI", "SII", "b4-b6", "a3", "a4", "a5", "L8")
memname_ras <- membership_ras
for (i in 1:length(community_ras)) { memname_ras[membership_ras==i] <- community_ras[i]}
i1 <- pdbs$resno[1, res3b[, 3]]; i2 <- pdbs$resno[1, res3b[, 4]]
ann1 <- memname_ras[i1]; ann2 <- memname_ras[i2]
res4b <- cbind(res3b[, 1:4], ann1, ann2, res3b[5:ncol(res3b)])

# add PC loadings into the table
lmt <- matrix(NA, ncol(pdbs$ali), ncol(pdbs$ali))
lmt[gaps.res$f.inds, gaps.res$f.inds] <- lm1
pc1 <- round(lmt[cbind(res4b$i, res4b$j)], 2)
lmt <- matrix(NA, ncol(pdbs$ali), ncol(pdbs$ali))
lmt[gaps.res$f.inds, gaps.res$f.inds] <- lm2
pc2 <- lmt[cbind(res4b$i, res4b$j)]
res4b <- cbind(res4b, pc1=pc1)
res4b <- cbind(res4b, pc2=pc2)

head(res4b, n=20)

# pymol.dccm(lm1, pdb=pdbs2pdb(pdbs, 1, rm.gaps=TRUE)[[1]], omit=0.24, step=0.01, file='pc_loading1.pml')

dm_ras <- dm.aa[gaps.res$f.inds, gaps.res$f.inds, ]
save(dm_ras,
     file="dm.RData")

save(pdbs.aa, gaps.res, gaps.pos, grps, res2b, res3b, res4b, 
     dm.ca, dm.aa, dm.aa.cut, dm.aa.masked, cm, xyz, 
     pc.ca, pc.aa, pc.aa.cut, pc.aa.masked, pc.cm, pc,
     file="eddm_ras.RData")


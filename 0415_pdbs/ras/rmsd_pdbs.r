## name: rmsd_pdbs.r
## date: 09/12/2016

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pca.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs_annotation.RData")

## rmsd between pdb structures
rmsd_pdbs <- rmsd(xyz_ras, a.inds=atom2xyz(which(!is.gap(pdbs_ras_new))), ncore=8)
rownames(rmsd_pdbs) <- pdbs_ras_new$id; colnames(rmsd_pdbs) <- pdbs_ras_new$id

## mean rmsd within/between groups
inds <- list("GTP" = which(tbl_ras[,"ligandColors"] == "red"),
             "GDP" = which(tbl_ras[,"ligandColors"] == "green"),
             "STATE1" = which(tbl_ras[,"ligandColors"] == "orange"),
             "GEF" = which(tbl_ras[,"ligandColors"] == "purple"))

rmsd_group <- matrix(NA, nrow=4, ncol=4)
rownames(rmsd_group) <- c("GTP", "GDP", "STATE1", "GEF")
colnames(rmsd_group) <- c("GTP", "GDP", "STATE1", "GEF")
diag(rmsd_group) <- c(mean(rmsd_pdbs[inds[[1]],inds[[1]]][lower.tri(rmsd_pdbs[inds[[1]],inds[[1]]])]),
                      mean(rmsd_pdbs[inds[[2]],inds[[2]]][lower.tri(rmsd_pdbs[inds[[2]],inds[[2]]])]),
                      mean(rmsd_pdbs[inds[[3]],inds[[3]]][lower.tri(rmsd_pdbs[inds[[3]],inds[[3]]])]),
                      mean(rmsd_pdbs[inds[[4]],inds[[4]]][lower.tri(rmsd_pdbs[inds[[4]],inds[[4]]])]))
for(i in 1:(dim(rmsd_group)[1]-1)) {
  for(j in (i+1):dim(rmsd_group)[1]) {
    rmsd_group[i,j] <- mean(rmsd_pdbs[inds[[i]],inds[[j]]])
  }
}

round(rmsd_group, digits=2)

save(rmsd_pdbs, rmsd_group, 
     file="rmsd_pdbs.RData")



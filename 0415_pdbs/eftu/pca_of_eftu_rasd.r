## name: pca_of_eftu_rasd.r
## date: 05/31/2016

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pca_of_eftu.RData")

# previous core fitting is based on RasD; here I directly perform PCA using previous fitting
# RasD boundary: 1D8T_B residues 9-198, pdbs_eftu$ali positions 1-212
inds1 <- 1:212
inds2 <- 1:(212*3)

pca_eftu_rasd <- pca(pdbs_eftu$xyz[,inds2], rm.gaps=TRUE)
# pca_eftu <- pca(pdbs_eftu$xyz, rm.gaps=TRUE)
# pca_eftu <- pca.xyz(xyz[, gaps.pos$f.inds])

plot.pca(pca_eftu_rasd, col=as.character(eftu_annotate[,"colorLig"]))
dev.copy2pdf(file="figures/pca_eftu_rasd_ligs.pdf")
plot.pca(pca_eftu_rasd, col=as.character(eftu_annotate[,"colorSource"]))
dev.copy2pdf(file="figures/pca_eftu_rasd_source.pdf")

# PCA on all structures
gaps.res <- gap.inspect(pdbs_eftu$ali[,inds1])
gaps.pos <- gap.inspect(pdbs_eftu$xyz[,inds2])

# Write PC trajectory
a <- mktrj.pca(pca_eftu_rasd, pc=1, file="pca/eftu_rasd_pc1.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )
b <- mktrj.pca(pca_eftu_rasd, pc=2, file="pca/eftu_rasd_pc2.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )

save(pca_eftu_rasd, pdbs_eftu, pdbs_eftu_new, eftu_annotate,
     file="pca_of_eftu_rasd.RData")

# 4P3Y_A is an outlier - SI is deformed I think.
plot(cbind(pca_eftu_rasd$z[,1],pca_eftu_rasd$z[,2]), cex=1.0, col=as.character(eftu_annotate[,"colorLig"]),
  xlab="pc1", ylab="pc2", main="pca_eftu_rasd_pc12_ligs")
identify(pca_eftu_rasd$z[,1],pca_eftu_rasd$z[,2],labels=substr(basename(pdbs_eftu$id),1,6))
dev.copy2pdf(file="figures/pc12_eftu_rasd_ligs.pdf")







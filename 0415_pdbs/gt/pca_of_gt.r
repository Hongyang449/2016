## name: pca_of_gt.r
## date: 06/05/2016

attach(transducin)
pc.xray <- pca(pdbs$xyz, rm.gaps=TRUE)
gaps.res <- gap.inspect(pdbs$ali)
gaps.pos <- gap.inspect(pdbs$xyz)


# Write PC trajectory
a <- mktrj.pca(pc.xray, pc=1, file="pca/gt_xray_pc1.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )
b <- mktrj.pca(pc.xray, pc=2, file="pca/gt_xray_pc2.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )
c <- mktrj.pca(pc.xray, pc=3, file="pca/gt_xray_pc3.pdb",
               resno = pdbs$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs$ali[1, gaps.res$f.inds]) )


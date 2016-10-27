## name: pca_loading_rasd.r
## date: 05/30/2016

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pca_of_eftu_rasd.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/sse_eftu.RData")

ylim=c(0,0.4)

# The variances captured by PCs
pca_eftu_rasd$L/sum(pca_eftu_rasd$L)

layout(matrix(1:2,nrow=2))
plot.bio3d(pca_eftu_rasd$au[,1],xaxt="n",typ="h",ylim=ylim,sse=sse_eftu,sse.min.length=3,
  xlab="Residue Number",ylab="PC1(85.73%)")
plot.bio3d(pca_eftu_rasd$au[,2],xaxt="n",typ="h",ylim=ylim,sse=sse_eftu,sse.min.length=3,
  xlab="Residue Number",ylab="PC2(8.47%)")
mtext("eftu_rasd",outer=T,line=-3,cex=1.5)
dev.copy2pdf(file="figures/pca_loading_eftu_rasd_sse.pdf")



## name: pca_loading.r
## date: 05/30/2016

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pca.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/sse_eftu.RData")

ylim=c(0,0.2)

# The variances captured by PCs
pca_eftu$L/sum(pca_eftu$L)

layout(matrix(1:2,nrow=2))
plot.bio3d(pca_eftu$au[,1],xaxt="n",typ="h",ylim=ylim,sse=sse_eftu,sse.min.length=3,
  xlab="Residue Number",ylab="PC1(94.11%)")
plot.bio3d(pca_eftu$au[,2],xaxt="n",typ="h",ylim=ylim,sse=sse_eftu,sse.min.length=3,
  xlab="Residue Number",ylab="PC2(3.59%)")
mtext("eftu",outer=T,line=-3,cex=1.5)
dev.copy2pdf(file="figures/pca_loading_eftu_sse.pdf")



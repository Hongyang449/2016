## name: pca_loading.r
## date: 04/18/2016

load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/gt/sse_gt.RData")

source("/Users/hyangl/project/ras/results/2015/functions/add.community.rect.R")

attach(transducin)
pc.xray <- pca(pdbs$xyz, rm.gaps=TRUE)
resno <- pdbs$resno[1,!is.gap(pdbs)]
ylim=c(0,0.4)

membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])

# The variances captured by PCs
pc.xray$L/sum(pc.xray$L)

#####################
## community label ##
#####################

layout(matrix(1:2,nrow=2))
plot(pc.xray$au[,1],xaxt="n",typ="h",ylim=ylim,xlab="Residue Number",ylab="PC1(49.79%)")
# add community partition
add.community.rect(ylim=ylim, membership=membership_gt)
num_community <- length(unique(membership_gt))
for (i in 1:num_community) {
  pos <- bounds(which(membership_gt == i))[,"start"]
  axis(1, at=pos, labels=resno[pos])
}
plot(pc.xray$au[,2],xaxt="n",typ="h",ylim=ylim,xlab="Residue Number",ylab="PC2(15.60%)")
# add community partition
add.community.rect(ylim=ylim, membership=membership_gt)
num_community <- length(unique(membership_gt))
for (i in 1:num_community) {
  pos <- bounds(which(membership_gt == i))[,"start"]
  axis(1, at=pos, labels=resno[pos])
}
mtext("gt",outer=T,line=-3,cex=1.5)
dev.copy2pdf(file="figures/pca_loading_gt_community.pdf")

###############
## sse label ##
###############

layout(matrix(1:2,nrow=2))
plot.bio3d(pc.xray$au[,1],xaxt="n",typ="h",ylim=ylim,sse=sse_gt,sse.min.length=3,
  xlab="Residue Number",ylab="PC1(51.23%)")
plot.bio3d(pc.xray$au[,2],xaxt="n",typ="h",ylim=ylim,sse=sse_gt,sse.min.length=3,
  xlab="Residue Number",ylab="PC2(18.80%)")
mtext("gt",outer=T,line=-3,cex=1.5)
dev.copy2pdf(file="figures/pca_loading_gt_sse.pdf")


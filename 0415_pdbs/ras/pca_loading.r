## name: pca_loading.r
## date: 04/18/2016

load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pca_plot.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/sse_ras.RData")

source("/Users/hyangl/project/ras/results/2015/functions/add.community.rect.R")

resno <- 1:166
ylim=c(0,0.4)

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

# The variances captured by PCs
pc.xray$L/sum(pc.xray$L)

#####################
## community label ##
#####################

layout(matrix(1:2,nrow=2))
plot(pc.xray$au[,1],xaxt="n",typ="h",ylim=ylim,xlab="Residue Number",ylab="PC1(51.23%)")
# add community partition
add.community.rect(ylim=ylim, membership=membership_ras)
num_community <- length(unique(membership_ras))
for (i in 1:num_community) {
  pos <- bounds(which(membership_ras == i))[,"start"]
  axis(1, at=pos, labels=resno[pos])
}
plot(pc.xray$au[,2],xaxt="n",typ="h",ylim=ylim,xlab="Residue Number",ylab="PC2(18.80%)")
# add community partition
add.community.rect(ylim=ylim, membership=membership_ras)
num_community <- length(unique(membership_ras))
for (i in 1:num_community) {
  pos <- bounds(which(membership_ras == i))[,"start"]
  axis(1, at=pos, labels=resno[pos])
}
mtext("ras",outer=T,line=-3,cex=1.5)
dev.copy2pdf(file="figures/pca_loading_ras_community.pdf")

###############
## sse label ##
###############

layout(matrix(1:2,nrow=2))
plot.bio3d(pc.xray$au[,1],xaxt="n",typ="h",ylim=ylim,sse=sse_ras,sse.min.length=3,
  xlab="Residue Number",ylab="PC1(51.23%)")
plot.bio3d(pc.xray$au[,2],xaxt="n",typ="h",ylim=ylim,sse=sse_ras,sse.min.length=3,
  xlab="Residue Number",ylab="PC2(18.80%)")
mtext("ras",outer=T,line=-3,cex=1.5)
dev.copy2pdf(file="figures/pca_loading_ras_sse.pdf")


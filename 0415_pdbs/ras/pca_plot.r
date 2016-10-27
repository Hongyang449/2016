## name: pca_plot.r
## date: 09/12/2016

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pca.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs_annotation.RData")

# 4 colors
cols <- tbl_ras[,"ligandColors"]
cols[cols=="green"] <- "green3"
table(cols)
# green3 orange purple    red
#     32      5     18     88
# isoform
cols_iso <- tbl_ras[,"isoform"]
cols_iso[tbl_ras[,"isoform"]=="H"] <- "blue"
cols_iso[tbl_ras[,"isoform"]=="K"] <- "green3"
# mutation
cols_mut <- tbl_ras[,"mutation"]
cols_mut[tbl_ras[,"mutation"]==""] <- "red"
cols_mut[tbl_ras[,"mutation"]!=""] <- "blue"
# extraLigand
cols_lig <- tbl_ras[,"extraLigand"]
cols_lig[tbl_ras[,"extraLigand"]=="F"] <- "red"
cols_lig[tbl_ras[,"extraLigand"]=="T"] <- "blue"
# hclust 
hc <- hclust(dist(pca_ras$z[, 1:3]))
# plot(hc, labels=pdbs_ras_new$id, main = "hclust_pc13_ras")
grps <- cutree(hc, k=4)
cols_h4 <- rep(NA, dim(tbl_ras)[1])
for(i in 1:length(cols_h4)) {
  if(grps[i] == 1) cols_h4[i] <- c("red")
  if(grps[i] == 2) cols_h4[i] <- c("purple")
  if(grps[i] == 3) cols_h4[i] <- c("green3")
  if(grps[i] == 4) cols_h4[i] <- c("orange")
}

plot(pca_ras, pc.axes=1:2, col=cols, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("GTP","GDP","GEF","STATE1"), 
  col=c("red","green3","purple","orange"), pch=1, cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_4_colors.pdf")

plot(pca_ras, pc.axes=1:2, col=cols_iso, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("HRAS","KRAS"), 
  col=c("blue","green3"), pch=1, cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_isoform.pdf")

plot(pca_ras, pc.axes=1:2, col=cols_mut, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("WT","MUTATION"), 
  col=c("red","blue"), pch=1, cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_mutation.pdf")

plot(pca_ras, pc.axes=1:2, col=cols_lig, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("NONE","EXTRA-LIGAND"), 
  col=c("red","blue"), pch=1, cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_extralig.pdf")

plot(pca_ras, pc.axes=1:2, col=cols_h4, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("GROUP1","GROUP2","GROUP3","GROUP4"),
  col=c("red","green3","purple","orange"), pch=1, cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_hcluster4.pdf")

layout(matrix(1:4, nrow=2))
plot(pca_ras, pc.axes=1:2, col=cols, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("GTP","GDP","GEF","STATE1"), 
  col=c("red","green3","purple","orange"), pch=1, cex=0.7)
plot(pca_ras, pc.axes=1:2, col=cols_iso, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("HRAS","KRAS"), 
  col=c("blue","green3"), pch=1, cex=0.7)
plot(pca_ras, pc.axes=1:2, col=cols_mut, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("WT","MUTATION"), 
  col=c("red","blue"), pch=1, cex=0.7)
plot(pca_ras, pc.axes=1:2, col=cols_lig, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("NONE","EXTRA-LIGAND"), 
  col=c("red","blue"), pch=1, cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_all.pdf")

## pc12 with rmsd

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/rmsd_pdbs.RData")

# coordinates of center for each group
pc1 <- c(mean(pca_ras$z[cols=="red",1]),mean(pca_ras$z[cols=="green3",1]),
         mean(pca_ras$z[cols=="orange",1]),mean(pca_ras$z[cols=="purple",1]))
pc2 <- c(mean(pca_ras$z[cols=="red",2]),mean(pca_ras$z[cols=="green3",2]),
         mean(pca_ras$z[cols=="orange",2]),mean(pca_ras$z[cols=="purple",2]))
# center between groups
mat1 <- matrix(NA, nrow=4, ncol=4); mat2 <- mat1
for (i in 1:(dim(mat1)[1]-1)) {
  for (j in (i+1):dim(mat1)[1]) {
    mat1[i,j] <- mean(c(pc1[i],pc1[j]))
    mat2[i,j] <- mean(c(pc2[i],pc2[j]))
  }
}
diag(mat1) <- pc1; diag(mat2) <- pc2

plot(pca_ras, pc.axes=1:2, col=cols, pch=15, xlim=c(-10,50), ylim=c(-10,25))
legend(x = "topright", legend = c("GTP","GDP","GEF","STATE1"),
  col=c("red","green3","purple","orange"), pch=1, cex=0.7)
# draw lines between group centers
lines(x=pc1[c(1,3,2,4,1,2,3,4)], y=pc2[c(1,3,2,4,1,2,3,4)], col="gray90", lty=2)
# add intra and inter groups rmsd
text(mat1[upper.tri(mat1, diag=T)], mat2[upper.tri(mat2, diag=T)],
  labels=round(rmsd_group, digits=2)[upper.tri(rmsd_group, diag=T)], cex=0.7)
dev.copy2pdf(file="figures/pc12_ras_4_colors_with_rmsd.pdf")




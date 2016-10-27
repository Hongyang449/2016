## name: net_ras_cmap_noh.r
## date: 01/26/2016

## Here I color edges that either significant diff OR more than 2 fold diff

library(abind)

load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/cmap_pdbs.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs_annotation.RData")

source("/Users/hyangl/project/ras/results/2015/functions/get.community.cij.R")
source("/Users/hyangl/project/ras/results/2015/functions/community.cij.R")
source("/Users/hyangl/project/ras/results/2015/functions/remodel.cna.R")
source("/Users/hyangl/project/ras/results/2015/functions/get.signif.R")
source("/Users/hyangl/project/ras/results/2015/functions/pairwise.signif.R")
source("/Users/hyangl/project/ras/results/2015/functions/pairwise.magnitude.R")

# !! here p.cutoff = 0.01 is used !! (cij network we use 0.05 and 2-fold)
p.cutoff <- 0.01
cutoff_magnitude=1/2
ifabs = T

cols <- tbl_ras[,"ligandColors"]
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

# build dummy networks!
cij_gtp <- filter.dccm(cmap_ras_noh[["6"]][,,cols == "red"], cutoff.cij=0)
cij_gdp <- filter.dccm(cmap_ras_noh[["6"]][,,cols == "green"], cutoff.cij=0)
cna_gtp <- cna(cij_gtp, cutoff.cij=0)
cna_gdp <- cna(cij_gdp, cutoff.cij=0)

net_complete <- lapply(cmap_ras_noh, function(x) {
  community_cij_ras_gtp <- get.community.cij(cij=x[,,cols == "purple"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  community_cij_ras_gdp <- get.community.cij(cij=x[,,cols == "green"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  # calculate p-value and whether magnitude diff larger than cutoff
  p_community_cij_ras <- pairwise.signif(list(gtp=community_cij_ras_gtp$raw,
    gdp=community_cij_ras_gdp$raw), method="wilcox.test")
  # remodel networks
  nets_ras_gtp_vs_gdp_signif <- remodel.cna(list(gtp=cna_gtp,gdp=cna_gdp),
    member=membership_ras, method="sum", col.edge="significance", scut=4,
    normalize=FALSE, signif=p_community_cij_ras[[1]], p.cutoff=p.cutoff)
  return(list(nets_ras_gtp_vs_gdp_signif, community_cij_ras_gtp, community_cij_ras_gdp,
    p_community_cij_ras))
  })

net_ras_cmap_noh <- lapply(net_complete, function (x) { x[[1]] })

save(net_complete, net_ras_cmap_noh,
     file="net_ras_cmap_noh.RData")


##### plot!
source("/Users/hyangl/project/ras/results/2015/functions/plot.nets.R")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/net_ras_cmap_noh.RData")
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/network/layout_2d.RData")
layout_2d <- layout_2d[1:9,]

plot.nets(net_ras_cmap_noh, layout_2d=layout_2d, width=0.1)
mtext("ras_cmap_noh_3_7", outer=T, line=-5)
dev.copy2pdf(file="figures/nets_ras_cmap_noh_3_7.pdf")

## new plot

source("/Users/hyangl/project/ras/results/2015/functions/trim.cna.R")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/net_ras_cmap_noh.RData")
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/network/layout_2d.RData")
layout_2d <- layout_2d[1:9,]
cutoff = 0.001
width = 0.1

layout(matrix(1:length(net_ras_cmap_noh), nrow=1))
invisible(lapply(as.list(names(net_ras_cmap_noh)), function(x) {
  net_ras <- trim.cna(net_ras_cmap_noh[[x]], cutoff_community=cutoff,
    membership = net_ras_cmap_noh[[x]][[1]]$communities$membership)
  plot.cna(net_ras, layout=layout_2d, w=(E(net_ras$community.network)$weight) * width, 
    vertex.label=NA, edge.label=NA)
  mtext(x, outer=F, line=-30)
  }))
mtext("ras_cmap_noh_3_7", outer=T, line=-27)
dev.copy2pdf(file="figures/nets_ras_cmap_noh_3_7_trimmed.pdf")















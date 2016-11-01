## name: net_ras_cmap_noh.r
## date: 01/26/2016

## color edges by significant diff of wilcox.test
## use the same cutoff as dynamic networks - p.cutoff 0.05

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

p.cutoff <- 0.05
ifabs = T

cols <- tbl_ras[,"ligandColors"]
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

net_complete <- lapply(cmap_ras_noh, function(x) {
  community_cij_ras_gtp <- get.community.cij(cij=x[,,cols == "red"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  community_cij_ras_gdp <- get.community.cij(cij=x[,,cols == "green"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  # calculate p-value and whether magnitude diff larger than cutoff
  p_community_cij_ras <- pairwise.signif(list(gtp=community_cij_ras_gtp$raw,
    gdp=community_cij_ras_gdp$raw), method="wilcox.test")
  # calculate average cij from different structure cmap
  cij_gtp <- filter.dccm(x[,,cols == "red"], cutoff.cij=0)
  cij_gdp <- filter.dccm(x[,,cols == "green"], cutoff.cij=0)  
  # build networks
  cna_gtp <- cna(cij_gtp, cutoff.cij=0)
  cna_gdp <- cna(cij_gdp, cutoff.cij=0)
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
load("/Users/hyangl/project/ras/results/2016/info/layout_2d.RData")

plot.nets(net_ras_cmap_noh, layout_2d=layout_ras, width=0.2)
mtext("ras_cmap_noh_3_6", outer=T, line=-5)
dev.copy2pdf(file="figures/nets_ras_cmap_noh_3_6.pdf")

## new plot

source("/Users/hyangl/project/ras/results/2015/functions/trim.cna.R")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/net_ras_cmap_noh.RData")
load("/Users/hyangl/project/ras/results/2016/info/layout_2d.RData")
cutoff = 0.001
width = 0.2

layout(matrix(1:length(net_ras_cmap_noh), nrow=1))
invisible(lapply(as.list(names(net_ras_cmap_noh)), function(x) {
  net_ras <- trim.cna(net_ras_cmap_noh[[x]], cutoff_community=cutoff,
    membership = net_ras_cmap_noh[[x]][[1]]$communities$membership)
  plot.cna(net_ras, layout=layout_ras, w=(E(net_ras$community.network)$weight) * width, 
    vertex.label=NA, edge.label=NA)
  mtext(x, outer=F, line=-30)
  }))
mtext("ras_cmap_noh_3_6", outer=T, line=-27)
dev.copy2pdf(file="figures/nets_ras_cmap_noh_3_6_trimmed.pdf")















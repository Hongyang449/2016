## name: net_dm.r
## date: 11/04/2016

## Here I want to build networks using distance matrix; the cij is defined as this:
## 1. 1/dm (for raw dm)
## 2. 1/masked_dm - 1/maksed_upperbound_distance (e.g. 8 or 10 here; for masked dm)


load("/Users/hyangl/project/ras/results/2016/1027_dm/ras/dm.RData")
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs_annotation.RData")

source("/Users/hyangl/bio3d/new_funs/eddm_funs.R")
source("/Users/hyangl/project/ras/results/2015/functions/get.community.cij.R")
source("/Users/hyangl/project/ras/results/2015/functions/community.cij.R")
source("/Users/hyangl/project/ras/results/2015/functions/remodel.cna.R")
source("/Users/hyangl/project/ras/results/2015/functions/get.signif.R")
source("/Users/hyangl/project/ras/results/2015/functions/pairwise.signif.R")

p.cutoff <- 0.05
ifabs = T

cols <- tbl_ras[,"ligandColors"]
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])

## raw dm
net_ras_raw_0.05 <- lapply(as.list(1), function(x) {
  cij <- 1/dm_ras
  community_cij_ras_gtp <- get.community.cij(cij=cij[,,cols == "red"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  community_cij_ras_gdp <- get.community.cij(cij=cij[,,cols == "green"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  # calculate p-value and whether magnitude diff larger than cutoff
  p_community_cij_ras <- pairwise.signif(list(gtp=community_cij_ras_gtp$raw,
    gdp=community_cij_ras_gdp$raw), method="wilcox.test")
  # calculate average cij from different structure cmap
  cij_gtp <- filter.dccm(cij[,,cols == "red"], cutoff.cij=0)
  cij_gdp <- filter.dccm(cij[,,cols == "green"], cutoff.cij=0)
  # build dummy networks
  cna_gtp <- cna(cij_gtp, cutoff.cij=0.4)
  cna_gdp <- cna(cij_gdp, cutoff.cij=0.4)
  # minus log cij
  cij_gtp[cij_gtp >= 0.9999] = 0.9999; cij_gtp <- -log(cij_gtp)
  cij_gtp[is.infinite(cij_gtp)] = 0; cij_gtp[lower.tri(cij_gtp)] <- t(cij_gtp)[lower.tri(cij_gtp)]
  cij_gdp[cij_gdp >= 0.9999] = 0.9999; cij_gdp <- -log(cij_gdp)
  cij_gdp[is.infinite(cij_gdp)] = 0; cij_gdp[lower.tri(cij_gdp)] <- t(cij_gdp)[lower.tri(cij_gdp)]
  # assign true cij (it's too slow to build cna using true cij)
  cna_gtp$cij <- cij_gtp
  cna_gdp$cij <- cij_gdp
  # remodel networks
  nets_ras_gtp_vs_gdp_signif <- remodel.cna(list(gtp=cna_gtp,gdp=cna_gdp),
    member=membership_ras, method="sum", col.edge="significance", scut=4,
    normalize=FALSE, signif=p_community_cij_ras[[1]], p.cutoff=p.cutoff)
  return(list(nets_ras_gtp_vs_gdp_signif, community_cij_ras_gtp, community_cij_ras_gdp,
    p_community_cij_ras))
  })
net_ras_dm_0.05 <- lapply(net_ras_raw_0.05, function (x) { x[[1]] })

## mask x - 8 ; larger than 8 is treated as 8; p=0.05
net_ras_8_0.05 <- lapply(as.list(seq(3,6,by=0.5)), function(x) {
  cij <- 1/mask.dm(dm_ras, x, 8) - 1/mask.dm(8, x, 8)
  community_cij_ras_gtp <- get.community.cij(cij=cij[,,cols == "red"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  community_cij_ras_gdp <- get.community.cij(cij=cij[,,cols == "green"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  # calculate p-value and whether magnitude diff larger than cutoff
  p_community_cij_ras <- pairwise.signif(list(gtp=community_cij_ras_gtp$raw,
    gdp=community_cij_ras_gdp$raw), method="wilcox.test")
  # calculate average cij from different structure cmap
  cij_gtp <- filter.dccm(cij[,,cols == "red"], cutoff.cij=0)
  cij_gdp <- filter.dccm(cij[,,cols == "green"], cutoff.cij=0)
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
net_ras_masked_dm_8_0.05 <- lapply(net_ras_8_0.05, function (x) { x[[1]] })

## mask x - 10 ; larger than 10 is treated as 10; p=0.05
net_ras_10_0.05 <- lapply(as.list(seq(3,6,by=0.5)), function(x) {
  cij <- 1/mask.dm(dm_ras, x, 10) - 1/mask.dm(10, x, 10)
  community_cij_ras_gtp <- get.community.cij(cij=cij[,,cols == "red"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  community_cij_ras_gdp <- get.community.cij(cij=cij[,,cols == "green"],
    member=membership_ras, cutoff.cij=0, abs=ifabs)
  # calculate p-value and whether magnitude diff larger than cutoff
  p_community_cij_ras <- pairwise.signif(list(gtp=community_cij_ras_gtp$raw,
    gdp=community_cij_ras_gdp$raw), method="wilcox.test")
  # calculate average cij from different structure cmap
  cij_gtp <- filter.dccm(cij[,,cols == "red"], cutoff.cij=0)
  cij_gdp <- filter.dccm(cij[,,cols == "green"], cutoff.cij=0)
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
net_ras_masked_dm_10_0.05 <- lapply(net_ras_10_0.05, function (x) { x[[1]] })

save(dm_ras,
     net_ras_dm_0.05, net_ras_masked_dm_8_0.05, net_ras_masked_dm_10_0.05,
     net_ras_raw_0.05, net_ras_8_0.05, net_ras_10_0.05,
     file="net_dm.RData")

## plot
source("/Users/hyangl/project/ras/results/2015/functions/plot.nets.R")
load("/Users/hyangl/project/ras/results/2016/1027_dm/ras/net_dm.RData")
load("/Users/hyangl/project/ras/results/2016/info/layout_2d.RData")

plot.nets(net_ras_dm_0.05, layout_2d=layout_ras, width=0.1)
mtext("ras_raw_dm (p=0.05)", outer=T, line=-3)
dev.copy2pdf(file="figures/nets_ras_raw_dm_0.05.pdf")

plot.nets(net_ras_masked_dm_8_0.05, layout_2d=layout_ras, width=1)
mtext("ras_masked_dm 3/6 to 8 (p=0.05)", outer=T, line=-5)
dev.copy2pdf(file="figures/nets_ras_masked_dm_8_0.05.pdf")

plot.nets(net_ras_masked_dm_10_0.05, layout_2d=layout_ras, width=1)
mtext("ras_masked_dm 3/6 to 10 (p=0.05)", outer=T, line=-5)
dev.copy2pdf(file="figures/nets_ras_masked_dm_10_0.05.pdf")








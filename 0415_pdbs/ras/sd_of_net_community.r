## name: sd_of_net_community.r
## date: 10/31/2016

## Here I want to examine the sd (variance) of community edges within and between states (gtp vs gdp)

library(abind)

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/net_ras_cmap_noh.RData")

## structure of net_complete:
## list(nets_ras_gtp_vs_gdp_signif, community_cij_ras_gtp, community_cij_ras_gdp, p_community_cij_ras)

## noh 4.5; calculate the sum of sd of all community edges
sum(apply(net_complete[['4.5']][[2]]$raw, 1:2, sd))
# [1] 40.88062 # within gtp
sum(apply(net_complete[['4.5']][[3]]$raw, 1:2, sd))
# [1] 41.30724 # within gdp
sum(apply(abind(net_complete[['4.5']][[2]]$raw,net_complete[['4.5']][[3]]$raw,along=3), 1:2, sd))
# [1] 50.63161 # combine gtp and gdp - much larger

## noh 6; calculate the sum of sd of all community edges
sum(apply(net_complete[['6']][[2]]$raw, 1:2, sd))
# [1] 65.65264 # within gtp
sum(apply(net_complete[['6']][[3]]$raw, 1:2, sd))
# [1] 68.9406  # within gdp
sum(apply(abind(net_complete[['6']][[2]]$raw,net_complete[['6']][[3]]$raw,along=3), 1:2, sd))
# [1] 81.83473 # combine gtp and gdp - much larger

# plot

load("/Users/hyangl/project/ras/results/2016/info/layout_2d.RData")

layout(matrix(1:2, nrow=1))

sd_gtp <- apply(net_complete[['6']][[2]]$raw, 1:2, sd)
sd_gtp[upper.tri(sd_gtp)] <- 0; weight <- sd_gtp[sd_gtp != 0]
plot.cna(net_ras_cmap_noh[["6"]]$gtp, layout=layout_ras, w=weight,
  vertex.label=NA, edge.label=NA, edge.color="gray")
mtext("gtp_noh_6", outer=F, line=-5)

sd_gdp <- apply(net_complete[['6']][[3]]$raw, 1:2, sd)
sd_gdp[upper.tri(sd_gdp)] <- 0; weight <- sd_gdp[sd_gdp != 0]
plot.cna(net_ras_cmap_noh[["6"]]$gdp, layout=layout_ras, w=weight,
  vertex.label=NA, edge.label=NA, edge.color="gray")
mtext("gdp_noh_6", outer=F, line=-5)

mtext("sd_of_cmap_noh_6", outer=T, line=-3)
dev.copy2pdf(file="figures/sd_ras_cmap_noh_6.pdf")



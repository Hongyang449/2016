#+ name, include=FALSE
## name: eddm_sca_pairs.r
## date: 11/11/2016 

#'---
#'title: "Ensemble Difference Distance Matrices(eddm) analysis of Ras Crystallographic Structures and Statistical Coupling Analysis(sca) of G protein family"
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(fig.path='figures/', dev='png', dev.args=list(type="cairo"), dpi=300)
library(bio3d)
library(png)

#+ preamble, include=FALSE, eval=FALSE
library(rmarkdown)

#+ load_pdbs
# load distance matrix, sca matrix and residue number/name/community of Ras
load("/Users/hyangl/project/ras/results/2016/1027_dm/ras/dmaa.RData")
load("/Users/hyangl/project/ras/results/2016/1027_dm/ras/eddm.RData")
load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/ras/sca_pairs.RData")
load("/Users/hyangl/project/ras/results/2016/info/res_ras.RData")

#'
#' ## 1. EDDM Analysis of Ras Crystallographic GTP-GDP State Structures 
#'
#' Residue pairs are listed if the difference of their mean masked distances between states 
#' are larger than 0.5 angstrom.
#'

#+ fig.cap='The distribution of the mean distances between all residue pairs'
dmaa_ras_gtp = apply(dmaa_ras166[,,grps==1],1:2,mean)
dmaa_ras_gdp = apply(dmaa_ras166[,,grps==3],1:2,mean)
layout(matrix(1:2, 1, 2))
hist(dmaa_ras_gtp, xlab="distance (angstrom)", main="the GTP-bound state")
hist(dmaa_ras_gdp, xlab="distance (angstrom)", main="the GDP-bound state")

#+ summary of eddm beta0.5
# The differences between GTP and GDP state
head(subseteddm_ras_gxp_beta0.5)
# summary of intra and inter community pairs
tbl_ras_gxp_beta0.5
# number of intra-community pairs
sum(diag(tbl_ras_gxp_beta0.5))
# number of inter-community pairs
sum(tbl_ras_gxp_beta0.5[upper.tri(tbl_ras_gxp_beta0.5)])

#+ SII - a3
# residue pairs identified between SII and a3
subseteddm_ras_gxp_beta0.5[subseteddm_ras_gxp_beta0.5[,"community_i"]=="SII" &
  subseteddm_ras_gxp_beta0.5[,"community_j"]=="a3",-c(3:6)]

#+ a1/b1-b3 - a5
# residue pairs identified between a1/b1-b3 and a5
subseteddm_ras_gxp_beta0.5[subseteddm_ras_gxp_beta0.5[,"community_i"]=="a1/b1-b3" &
  subseteddm_ras_gxp_beta0.5[,"community_j"]=="a5", -(3:6)]

#'
#' ## 2. Comparison of EDDM and SCA results 
#'
#' The original SCA scores of all residue pairs are used.
#' The top 100 residue pairs from EDDM and SCA are compared. 
#'

#+ summary of eddm 100
# summary of eddm top 100
tbl_eddm_ras_gxp_100
# number of intra-community pairs
sum(diag(tbl_sca_ras_100))
# number of inter-community pairs
sum(tbl_sca_ras_100[upper.tri(tbl_sca_ras_100)])

#+ summary of sca 100
# summary of sca top 100
tbl_sca_ras_100
# number of intra-community pairs
sum(diag(tbl_sca_ras_100))
# number of inter-community pairs
sum(tbl_sca_ras_100[upper.tri(tbl_sca_ras_100)])

#+ common pairs found by sca and eddm
# total number of common pairs out of top 100
sum(pairs_sca %in% pairs_eddm)
sca_ras[(1:100)[pairs_sca %in% pairs_eddm],]
subseteddm_ras_gxp_beta0.5[(1:100)[pairs_eddm %in% pairs_sca],-c(3:4)]

#+ smaller data.frame summarizing common pairs
# the 10 common pairs
inds <- match(pairs_eddm, pairs_sca)[!is.na(match(pairs_eddm, pairs_sca))]
sca_ras[(1:100)[inds],]
common_pairs <- cbind(subseteddm_ras_gxp_beta0.5[(1:100)[pairs_eddm %in% pairs_sca],c(1:2,5:9)],
  sca=sca_ras[(1:100)[inds],"sca"])

#+ fisher.test
# fisher's exact test of the assocation between eddm and sca ranking.
fisher.test(matrix(c(10,90,90,l*(l-1)/2-100), 2,2))

#+ close, include=TRUE, eval=FALSE
library(rmarkdown)
render("eddm_sca_pairs.r", "pdf_document", clean=TRUE)

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)



#+ name, include=FALSE
## name: eddm_sca_pairs_gt.r
## date: 11/11/2016 

#'---
#'title: "Ensemble Difference Distance Matrices(eddm) analysis of Gat Crystallographic Structures and Statistical Coupling Analysis(sca) of Gat protein family (pfam PF00503)"
#'---

#+ setup, include=FALSE
knitr::opts_chunk$set(fig.path='figures/', dev='png', dev.args=list(type="cairo"), dpi=300)
library(bio3d)
library(png)

#+ preamble, include=FALSE, eval=FALSE
library(rmarkdown)

#+ load_pdbs
# load distance matrix, sca matrix and residue number/name/community of Gat
load("/Users/hyangl/project/ras/results/2016/1027_dm/gt/dmaa.RData")
load("/Users/hyangl/project/ras/results/2016/1027_dm/gt/eddm.RData")
load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/gt/sca_pairs.RData")
load("/Users/hyangl/project/ras/results/2016/info/res_gt.RData")

#'
#' ## 1. EDDM Analysis of Gat Crystallographic GTP-GDP State Structures 
#'
#' Residue pairs are listed if the difference of their mean masked distances between states 
#' are larger than 0.5 angstrom.
#'

#+ basic statistics of all-atom distance matrix
# general statistics of all-atom distance matrix of gt structures in GTP/GDP states
attach(transducin)
grps <- annotation[, 'state3']
grps[grps=='GTP'] <- 1; grps[grps=='GDP'] <- 2; grps[grps=='GDI'] <- 3;
dmaa_gt_gtp = apply(dmaa_gt305[,,grps==1],1:2,mean)
dmaa_gt_gdp = apply(dmaa_gt305[,,grps==2],1:2,mean)
# mean distance of all residue pairs
mean(dmaa_gt_gtp[upper.tri(dmaa_gt_gtp)]); mean(dmaa_gt_gdp[upper.tri(dmaa_gt_gdp)]);
# sd of all the distances
sd(dmaa_gt_gtp[upper.tri(dmaa_gt_gtp)]); sd(dmaa_gt_gdp[upper.tri(dmaa_gt_gdp)]);

#+ distribution of residue-residue distances, fig.cap='The distribution of the mean distances between all residue pairs'
layout(matrix(1:2, 1, 2))
hist(dmaa_gt_gtp, xlab="distance (angstrom)", main="the GTP-bound state")
hist(dmaa_gt_gdp, xlab="distance (angstrom)", main="the GDP-bound state")

#+ summary of eddm beta0.5
# The differences between GTP and GDP state
head(subseteddm_gt_gxp_beta0.5)
# summary of intra and inter community pairs
tbl_gt_gxp_beta0.5
# number of intra- and inter-community pairs
sum(diag(tbl_gt_gxp_beta0.5)); sum(tbl_gt_gxp_beta0.5[upper.tri(tbl_gt_gxp_beta0.5)])

#+ SII - a3
# residue pairs identified between SII and a3
subseteddm_gt_gxp_beta0.5[subseteddm_gt_gxp_beta0.5[,"community_i"]=="SII" &
  subseteddm_gt_gxp_beta0.5[,"community_j"]=="a3",-c(3:6)]

#+ a1/b1-b3 - a5
# residue pairs identified between a1/b1-b3 and a5
subseteddm_gt_gxp_beta0.5[subseteddm_gt_gxp_beta0.5[,"community_i"]=="a1/b1-b3" &
  subseteddm_gt_gxp_beta0.5[,"community_j"]=="a5", -(3:6)]

#+ SII - a3 from sca top 100
# residue pairs between SII and a3 from top 100 sca
sca_gt_100 <- sca_gt[1:100,]
sca_gt_100[sca_gt_100[,"community_i"]=="SII" & sca_gt_100[,"community_j"]=="a3",]
# residue pairs between a1/b1-b3 - a5
sca_gt_100[sca_gt_100[,"community_i"]=="a1/b1-b3" & sca_gt_100[,"community_j"]=="a5",]

#'
#' ## 2. Comparison of EDDM and SCA results 
#'
#' EDDM and SCA (PF00503) return residue pairs of the entire protein.
#' The top 100 residue pairs from EDDM and SCA are compared. 
#'

#+ summary of eddm 100
# summary of eddm top 100
tbl_eddm_gt_gxp_100
# number of intra- and inter-community pairs
sum(diag(tbl_eddm_gt_gxp_100)); sum(tbl_eddm_gt_gxp_100[upper.tri(tbl_eddm_gt_gxp_100)])

#+ summary of sca 100
# summary of sca top 100
tbl_sca_gt_100
# number of intra- and inter-community pairs
sum(diag(tbl_sca_gt_100)); sum(tbl_sca_gt_100[upper.tri(tbl_sca_gt_100)])

#+ common pairs found by sca and eddm
# total number of common pairs out of top 100
sum(pairs_sca %in% pairs_eddm)
sca_gt[(1:100)[pairs_sca %in% pairs_eddm],]
subseteddm_gt_gxp_beta0.5[(1:100)[pairs_eddm %in% pairs_sca],-c(3:4)]

#+ smaller data.frame summarizing common pairs
# the common pairs
common_pairs

#+ fisher.test
# fisher's exact test of the assocation between eddm and sca ranking. 
fisher.test(matrix(c( sum(pairs_sca %in% pairs_eddm), 100-sum(pairs_sca %in% pairs_eddm),
  100-sum(pairs_sca %in% pairs_eddm), l*(l-1)/2-100 ), 2,2))

#+ close, include=TRUE, eval=FALSE
library(rmarkdown)
render("eddm_sca_pairs_gt.r", "pdf_document", clean=TRUE)

#'
#' ## Information About the Current Bio3D Session
#'
print(sessionInfo(), FALSE)



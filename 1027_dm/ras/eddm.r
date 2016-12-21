## name: eddm.r
## date: 11/08/2016

## Here I perform eddm analysis of all-atom distance matrices of ras pdbs

load("/Users/hyangl/project/ras/results/2016/1027_dm/ras/dmaa.RData")
load("/Users/hyangl/project/ras/results/2016/info/res_ras.RData")

eddm_ras <- eddm(pdbs_ras_aa, dm=dmaa_ras, grps=grps)
# aa321
eddm_ras[,c("a","b")] <- t(apply(as.matrix(eddm_ras[,c("a","b")]),1,function(x) {
  c(gsub("^.{0,3}",aa321(substr(x[1],1,3)),x[1]),gsub("^.{0,3}",aa321(substr(x[2],1,3)),x[2]))
  }))

## difference in the mean masked distance >= 1.0 (beta=1.0)
# grps==1/3 are gtp/gdp states
subseteddm_ras_gxp_beta1.0 <- subset.eddm(eddm_ras, grps=c(1,3), alpha=0.05, beta=1.0)
# add sse annotations to the table
i1 <- pdbs_ras_aa$resno[1, subseteddm_ras_gxp_beta1.0 [, "i"]]
i2 <- pdbs_ras_aa$resno[1, subseteddm_ras_gxp_beta1.0 [, "j"]]
tmp = names(membership_ras); names(tmp) = resno_ras
community_i <- tmp[as.character(i1)]; community_j <- tmp[as.character(i2)]
subseteddm_ras_gxp_beta1.0  <- cbind(subseteddm_ras_gxp_beta1.0 [, 1:4], community_i, community_j,
  subseteddm_ras_gxp_beta1.0 [5:ncol(subseteddm_ras_gxp_beta1.0 )])

# summary the number of highlighted pairs between communities
tmp <- subseteddm_ras_gxp_beta1.0[, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(subseteddm_ras_gxp_beta1.0[, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_ras_gxp_beta1.0 <- matrix(0, length(community_ras), length(community_ras))
rownames(tbl_ras_gxp_beta1.0) <- community_ras; colnames(tbl_ras_gxp_beta1.0) <- community_ras;
tbl_ras_gxp_beta1.0[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_ras_gxp_beta1.0) <- diag(tbl_ras_gxp_beta1.0)/2
tbl_ras_gxp_beta1.0
#         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8
#a1/b1-b3        0  1  5   0     0  0  0  0  0
#PL              1  0  6   6     0  0  0  0  0
#SI              5  6  3  13     0  0  0  0  0
#SII             0  6 13  11     0  5  0  0  0
#b4-b6           0  0  0   0     1  0  0  0  0
#a3              0  0  0   5     0  1  0  0  0
#a4              0  0  0   0     0  0  0  0  0
#a5              0  0  0   0     0  0  0  0  0
#L8              0  0  0   0     0  0  0  0  0
sum(diag(tbl_ras_gxp_beta1.0))
#[1] 16
sum(tbl_ras_gxp_beta1.0[upper.tri(tbl_ras_gxp_beta1.0)])
#[1] 36

# SII - a3
subseteddm_ras_gxp_beta1.0[subseteddm_ras_gxp_beta1.0[,"community_i"]=="SII" & 
  subseteddm_ras_gxp_beta1.0[,"community_j"]=="a3",-c(3:6)]
#         a   b  d.1  d.3 dm.1_3 dz.1_3 md.1 md.3 mdm.1_3 mdz.1_3
#73-101 R68 Y96 3.80 7.37   3.56   2.84 3.71 6.28    2.57    2.93
#65-101 G60 Y96 3.92 6.74   2.82   4.09 3.90 6.07    2.17    3.59
#74-104 D69 Q99 3.75 6.46   2.71   2.43 3.71 5.81    2.09    2.02
#64-101 A59 Y96 4.89 7.89   3.01   3.25 4.79 6.44    1.65    2.39
#70-104 S65 Q99 5.19 8.15   2.96   2.23 4.99 6.46    1.46    1.52

## difference in the mean masked distance >= 0.5 (beta=0.5)
# grps==1/3 are gtp/gdp states
subseteddm_ras_gxp_beta0.5 <- subset.eddm(eddm_ras, grps=c(1,3), alpha=0.05, beta=0.5)
# add sse annotations to the table
i1 <- pdbs_ras_aa$resno[1, subseteddm_ras_gxp_beta0.5 [, "i"]]
i2 <- pdbs_ras_aa$resno[1, subseteddm_ras_gxp_beta0.5 [, "j"]]
tmp = names(membership_ras); names(tmp) = resno_ras
community_i <- tmp[as.character(i1)]; community_j <- tmp[as.character(i2)]
subseteddm_ras_gxp_beta0.5  <- cbind(subseteddm_ras_gxp_beta0.5 [, 1:4], community_i, community_j,
  subseteddm_ras_gxp_beta0.5 [5:ncol(subseteddm_ras_gxp_beta0.5 )])

# summary the number of highlighted pairs between communities
tmp <- subseteddm_ras_gxp_beta0.5[, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(subseteddm_ras_gxp_beta0.5[, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_ras_gxp_beta0.5 <- matrix(0, length(community_ras), length(community_ras))
rownames(tbl_ras_gxp_beta0.5) <- community_ras; colnames(tbl_ras_gxp_beta0.5) <- community_ras;
tbl_ras_gxp_beta0.5[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_ras_gxp_beta0.5) <- diag(tbl_ras_gxp_beta0.5)/2
tbl_ras_gxp_beta0.5
#            a1/b1-b3 PL SI SII b4-b6 a3 L8 a4 a5
#   a1/b1-b3        9  1  7  11     0  0  0  0  4
#   PL              1  0  9  15     0  0  0  0  0
#   SI              7  9 11  27     0  0  0  0  0
#   SII            11 15 27  27     0 13  0  0  0
#   b4-b6           0  0  0   0     2  1  2  0  0
#   a3              0  0  0  13     1  2  0  2  0
#   L8              0  0  0   0     2  0  0  0  0
#   a4              0  0  0   0     0  2  0  0  0
#   a5              4  0  0   0     0  0  0  0  0
sum(diag(tbl_ras_gxp_beta0.5))
#[1] 51 intra-community
sum(tbl_ras_gxp_beta0.5[upper.tri(tbl_ras_gxp_beta0.5)])
#[1] 92 inter-community

# Here we cannot say eddm highlights more inter-group pairs than intra-group pairs!
# In fact we find more intra pairs - fisher.test
sum(apply(table(membership_ras),1, function(x) { x*(x-1)/2 }))
#[1] 1838 total number of intra pairs
length(membership_ras) * (length(membership_ras)-1) /2
#[1] 13695 total number of intra and inter pairs

# SII and a3
subseteddm_ras_gxp_beta0.5[subseteddm_ras_gxp_beta0.5[,"community_i"]=="SII" &
  subseteddm_ras_gxp_beta0.5[,"community_j"]=="a3",-c(3:6)]
#         a    b   d.1   d.3 dm.1_3 dz.1_3 md.1 md.3 mdm.1_3 mdz.1_3
#73-101 R68  Y96  3.80  7.37   3.56   2.84 3.71 6.28    2.57    2.93
#65-101 G60  Y96  3.92  6.74   2.82   4.09 3.90 6.07    2.17    3.59
#74-104 D69  Q99  3.75  6.46   2.71   2.43 3.71 5.81    2.09    2.02
#64-101 A59  Y96  4.89  7.89   3.01   3.25 4.79 6.44    1.65    2.39
#70-104 S65  Q99  5.19  8.15   2.96   2.23 4.99 6.46    1.46    1.52
#74-107 D69 R102  5.08  6.97   1.90   1.09 4.81 5.79    0.98    0.70
#73-104 R68  Q99  3.45  4.67   1.22   1.29 3.40 4.36    0.96    1.37
#76-101 Y71  Y96  7.24  8.79   1.55   0.46 5.39 6.29    0.91    0.82
#77-101 M72  Y96  4.58  5.61   1.03   1.68 4.55 5.42    0.87    1.55
#66-101 Q61  Y96  5.49  7.12   1.63   1.78 5.34 6.07    0.74    0.94
#71-104 A66  Q99  6.71 10.34   3.63   2.66 5.96 6.55    0.59    0.81
#69-107 Y64 R102 11.57 10.34   1.23   0.45 6.43 5.88    0.55    1.23
#73-100 R68  Q95  6.54  8.02   1.48   0.83 5.80 6.30    0.50    0.53

# a1/b1-b3 - a5
subseteddm_ras_gxp_beta0.5[subseteddm_ras_gxp_beta0.5[,"community_i"]=="a1/b1-b3" &
  subseteddm_ras_gxp_beta0.5[,"community_j"]=="a5", -(3:6)]
#         a    b  d.1  d.3 dm.1_3 dz.1_3 md.1 md.3 mdm.1_3 mdz.1_3
#47-154 K42 R149 6.91 5.55   1.36   1.06 6.07 5.26    0.81    1.50
#55-169 T50 R164 6.62 5.70   0.93   1.17 6.08 5.46    0.62    1.68
#47-158 K42 E153 7.40 5.86   1.53   1.18 6.24 5.64    0.60    1.28
#56-169 C51 R164 4.71 4.14   0.57   0.63 4.65 4.13    0.52    0.66


save(eddm_ras, 
     subseteddm_ras_gxp_beta1.0, tbl_ras_gxp_beta1.0,
     subseteddm_ras_gxp_beta0.5, tbl_ras_gxp_beta0.5,
     file="eddm.RData")












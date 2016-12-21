## name: sca_pairs.r
## date: 11/11/2016

## Here I analyze the original sca matrix obtained from pySCA analysis and compare it with eddm pairs

load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/eftu/pfam/sca.RData")
load("/Users/hyangl/project/ras/results/2016/1027_dm/eftu/eddm.RData")
load("/Users/hyangl/project/ras/results/2016/info/res_eftu.RData")

source("/Users/hyangl/project/ras/results/2015/functions/position_vec2mat.R")

l = length(inds_eftu)
# dim 167*167; here sca only consider RasD
rownames(scaij_eftu) <- resnano_eftu[resno_eftu %in% inds_eftu]
colnames(scaij_eftu) <- resnano_eftu[resno_eftu %in% inds_eftu]
cij <- matrix(0, nrow=nrow(scaij_eftu), ncol=ncol(scaij_eftu)); cij[upper.tri(cij)] <- scaij_eftu[upper.tri(scaij_eftu)];
rownames(cij) <- resnano_eftu[resno_eftu %in% inds_eftu]
colnames(cij) <- resnano_eftu[resno_eftu %in% inds_eftu]

# full ranked sca pairs
inds <- position_vec2mat(order(cij, decreasing=T)[1:(l*(l-1)/2)], dim=l)
sca_eftu <- as.data.frame(t(apply(inds, 1, function(x) { c(rownames(cij)[x],round(cij[x[1],x[2]],2)) })))
# add community
tmp = names(membership_eftu); names(tmp) = resnano_eftu
community_i <- tmp[as.character(sca_eftu[,1])]; community_j <- tmp[as.character(sca_eftu[,2])]
sca_eftu <- cbind(sca_eftu[,1:2], community_i, community_j, sca_eftu[,3:ncol(sca_eftu)])
colnames(sca_eftu) <- c("a","b","community_i","community_j","sca")

# summary of inter and intra community pairs in top 100 pairs within RasD
tmp <- sca_eftu[1:100, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(sca_eftu[1:100, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_sca_eftu_100 <- matrix(0, length(community_eftu), length(community_eftu))
rownames(tbl_sca_eftu_100) <- community_eftu; colnames(tbl_sca_eftu_100) <- community_eftu; 
tbl_sca_eftu_100[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_sca_eftu_100) <- diag(tbl_sca_eftu_100)/2
tbl_sca_eftu_100
         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8 D2 D3
a1/b1-b3        3  4  7  10     0  1  0  0  0  0  0
PL              4  5  4  11     3  9  0  0  0  0  0
SI              7  4  9   8     0  3  0  0  0  0  0
SII            10 11  8  10     2  8  0  1  1  0  0
b4-b6           0  3  0   2     0  0  0  0  0  0  0
a3              1  9  3   8     0  0  0  0  1  0  0
a4              0  0  0   0     0  0  0  0  0  0  0
a5              0  0  0   1     0  0  0  0  0  0  0
L8              0  0  0   1     0  1  0  0  0  0  0
D2              0  0  0   0     0  0  0  0  0  0  0
D3              0  0  0   0     0  0  0  0  0  0  0
sum(diag(tbl_sca_eftu_100))
#[1] 27
sum(tbl_sca_eftu_100[upper.tri(tbl_sca_eftu_100)])
#[1] 73

# top 100 of eddm within RasD
inds_d2d3 <- (subseteddm_eftu_gxp_beta0.5[,"community_i"]=="D2" |
  subseteddm_eftu_gxp_beta0.5[,"community_i"]=="D3" |
  subseteddm_eftu_gxp_beta0.5[,"community_j"]=="D2" |
  subseteddm_eftu_gxp_beta0.5[,"community_j"]=="D3")
tmp <- subseteddm_eftu_gxp_beta0.5[which(!inds_d2d3)[1:100], c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(subseteddm_eftu_gxp_beta0.5[which(!inds_d2d3)[1:100], c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_eddm_eftu_gxp_100 <- matrix(0, length(community_eftu), length(community_eftu))
rownames(tbl_eddm_eftu_gxp_100) <- community_eftu; colnames(tbl_eddm_eftu_gxp_100) <- community_eftu;
tbl_eddm_eftu_gxp_100[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_eddm_eftu_gxp_100) <- diag(tbl_eddm_eftu_gxp_100)/2
tbl_eddm_eftu_gxp_100
         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8 D2 D3
a1/b1-b3        0  0  7   6     0  0  0  1  0  0  0
PL              0  0  1   4     0  0  0  0  0  0  0
SI              7  1 23  19     0  0  0  0  0  0  0
SII             6  4 19  21     0 10  0  0  0  0  0
b4-b6           0  0  0   0     0  0  0  0  0  0  0
a3              0  0  0  10     0  0  0  0  0  0  0
a4              0  0  0   0     0  0  0  0  0  0  0
a5              1  0  0   0     0  0  0  8  0  0  0
L8              0  0  0   0     0  0  0  0  0  0  0
D2              0  0  0   0     0  0  0  0  0  0  0
D3              0  0  0   0     0  0  0  0  0  0  0
sum(diag(tbl_eddm_eftu_gxp_100))
#[1] 52
sum(tbl_eddm_eftu_gxp_100[upper.tri(tbl_eddm_eftu_gxp_100)])
#[1] 48

# top 100 pairs
pairs_eddm <- apply(as.matrix(subseteddm_eftu_gxp_beta0.5[which(!inds_d2d3)[1:100],c("a","b")]),1,
  function(x) { paste0(x[1],"-",x[2]) }); names(pairs_eddm) <- 1:100
pairs_sca <- apply(as.matrix(sca_eftu[1:100,c("a","b")]),1,
  function(x) { paste0(x[1],"-",x[2]) }); 

# common pairs found by sca and eddm
sum(pairs_sca %in% pairs_eddm)
#[1] 8
sca_eftu[(1:100)[pairs_sca %in% pairs_eddm],]
subseteddm_eftu_gxp_beta0.5[which(!inds_d2d3)[1:100][pairs_eddm %in% pairs_sca],-c(3:4)]

# smaller data.frame summarizing common pairs
inds <- match(pairs_eddm, pairs_sca)[!is.na(match(pairs_eddm, pairs_sca))]
sca_eftu[(1:100)[inds],]
common_pairs <- cbind(subseteddm_eftu_gxp_beta0.5[which(!inds_d2d3)[1:100][pairs_eddm %in% pairs_sca],c(1:2,5:9)],
  sca=sca_eftu[(1:100)[inds],"sca"])

# fisher.test
fisher.test(matrix(c(8,92,92,l*(l-1)/2-100), 2,2))
#p-value = 5.67e-07

save(scaij_eftu, inds_eftu, pairs_eddm, pairs_sca, l,
     sca_eftu, tbl_sca_eftu_100, tbl_eddm_eftu_gxp_100, common_pairs, 
     file="sca_pairs.RData")



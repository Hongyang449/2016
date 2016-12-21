## name: sca_pairs.r
## date: 11/11/2016

## Here I analyze the original sca matrix obtained from pySCA analysis and compare it with eddm pairs

load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/gt/pfam/sca.RData")
load("/Users/hyangl/project/ras/results/2016/1027_dm/gt/eddm.RData")
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")

source("/Users/hyangl/project/ras/results/2015/functions/position_vec2mat.R")

l = length(inds_gt)
# dim 310*310
resno_gt <- as.numeric(ali["resno_gt",ali["resno_gt",]!="0"])
names(resno_gt) <- ali["1TND_A",ali["resno_gt",]!="0"]
resnano_gt <- paste0(names(resno_gt), resno_gt)

rownames(scaij_gt) <- resnano_gt[resno_gt %in% inds_gt]
colnames(scaij_gt) <- resnano_gt[resno_gt %in% inds_gt]
cij <- matrix(0, nrow=nrow(scaij_gt), ncol=ncol(scaij_gt)); cij[upper.tri(cij)] <- scaij_gt[upper.tri(scaij_gt)];
rownames(cij) <- resnano_gt[resno_gt %in% inds_gt]
colnames(cij) <- resnano_gt[resno_gt %in% inds_gt]

# full ranked sca pairs
inds <- position_vec2mat(order(cij, decreasing=T)[1:(l*(l-1)/2)], dim=l)
sca_gt <- as.data.frame(t(apply(inds, 1, function(x) { c(rownames(cij)[x],round(cij[x[1],x[2]],2)) })))
# add community
load("/Users/hyangl/project/ras/results/2016/info/res_gt.RData")
tmp = names(membership_gt); names(tmp) = resnano_gt
community_i <- tmp[as.character(sca_gt[,1])]; community_j <- tmp[as.character(sca_gt[,2])]
sca_gt <- cbind(sca_gt[,1:2], community_i, community_j, sca_gt[,3:ncol(sca_gt)])
colnames(sca_gt) <- c("a","b","community_i","community_j","sca")

# summary of inter and intra community pairs in top 100 pairs within RasD
tmp <- sca_gt[1:100, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(sca_gt[1:100, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_sca_gt_100 <- matrix(0, length(community_gt), length(community_gt))
rownames(tbl_sca_gt_100) <- community_gt; colnames(tbl_sca_gt_100) <- community_gt; 
tbl_sca_gt_100[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_sca_gt_100) <- diag(tbl_sca_gt_100)/2
tbl_sca_gt_100
         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8 HD1 HD2
a1/b1-b3        1  3  0   6     0  3  0  0  0   1   0
PL              3  2  0   4     0  0  0  0  0   0   0
SI              0  0  1   6     2  2  0  0  0   0   0
SII             6  4  6  18     8 17  0  0  3   5   0
b4-b6           0  0  2   8     3  4  0  1  2   1   0
a3              3  0  2  17     4  2  0  0  2   2   0
a4              0  0  0   0     0  0  0  0  0   0   0
a5              0  0  0   0     1  0  0  0  0   0   0
L8              0  0  0   3     2  2  0  0  0   1   0
HD1             1  0  0   5     1  2  0  0  1   0   0
HD2             0  0  0   0     0  0  0  0  0   0   0
sum(diag(tbl_sca_gt_100))
#[1] 27
sum(tbl_sca_gt_100[upper.tri(tbl_sca_gt_100)])
#[1] 73

# top 100 of eddm
tmp <- subseteddm_gt_gxp_beta0.5[1:100, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(subseteddm_gt_gxp_beta0.5[1:100, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_eddm_gt_gxp_100 <- matrix(0, length(community_gt), length(community_gt))
rownames(tbl_eddm_gt_gxp_100) <- community_gt; colnames(tbl_eddm_gt_gxp_100) <- community_gt;
tbl_eddm_gt_gxp_100[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_eddm_gt_gxp_100) <- diag(tbl_eddm_gt_gxp_100)/2
tbl_eddm_gt_gxp_100
         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8 HD1 HD2
a1/b1-b3        1  0  0  14     0  0  0  0  0   0   0
PL              0  0  1  14     0  1  0  0  0   0   0
SI              0  1  0  12     0  0  0  0  0   0   0
SII            14 14 12  21     1 21  0  0  0   0   0
b4-b6           0  0  0   1     0  0  0  0  0   0   0
a3              0  1  0  21     0  3  1  0  0   0   2
a4              0  0  0   0     0  1  1  0  0   0   0
a5              0  0  0   0     0  0  0  0  0   0   0
L8              0  0  0   0     0  0  0  0  2   0   0
HD1             0  0  0   0     0  0  0  0  0   0   3
HD2             0  0  0   0     0  2  0  0  0   3   2
sum(diag(tbl_eddm_gt_gxp_100))
#[1] 30
sum(tbl_eddm_gt_gxp_100[upper.tri(tbl_eddm_gt_gxp_100)])
#[1] 70

# top 100 pairs
pairs_eddm <- apply(as.matrix(subseteddm_gt_gxp_beta0.5[1:100,c("a","b")]),1,
  function(x) { paste0(x[1],"-",x[2]) }); names(pairs_eddm) <- 1:100
pairs_sca <- apply(as.matrix(sca_gt[1:100,c("a","b")]),1,
  function(x) { paste0(x[1],"-",x[2]) }); 

# common pairs found by sca and eddm
sum(pairs_sca %in% pairs_eddm)
#[1] 12
sca_gt[(1:100)[pairs_sca %in% pairs_eddm],]
subseteddm_gt_gxp_beta0.5[c(1:100)[pairs_eddm %in% pairs_sca],-c(3:4)]

# smaller data.frame summarizing common pairs
inds <- match(pairs_eddm, pairs_sca)[!is.na(match(pairs_eddm, pairs_sca))]
sca_gt[(1:100)[inds],]
common_pairs <- cbind(subseteddm_gt_gxp_beta0.5[c(1:100)[pairs_eddm %in% pairs_sca],c(1:2,5:9)],
  sca=sca_gt[(1:100)[inds],"sca"])

# fisher.test
fisher.test(matrix(c(12,88,88,l*(l-1)/2-100), 2,2))
#p-value = 2.2e-16

save(scaij_gt, inds_gt, pairs_eddm, pairs_sca, l,
     sca_gt, tbl_sca_gt_100, tbl_eddm_gt_gxp_100, common_pairs, 
     file="sca_pairs.RData")



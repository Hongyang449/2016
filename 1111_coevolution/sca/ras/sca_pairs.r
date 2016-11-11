## name: sca_pairs.r
## date: 11/11/2016

## Here I analyze the original sca matrix obtained from pySCA analysis and compare it with eddm pairs

load("/Users/hyangl/project/ras/results/2016/1111_coevolution/sca/ras/sca.RData")
load("/Users/hyangl/project/ras/results/2016/1027_dm/ras/eddm.RData")
load("/Users/hyangl/project/ras/results/2016/info/res_ras.RData")


inds_sca <- c('5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157', '158', '159', '160', '161', '162', '165')
l = length(inds_sca)
# dim 158*158
rownames(cijr) <- resnano_ras[resno_ras %in% inds_sca]
colnames(cijr) <- resnano_ras[resno_ras %in% inds_sca]
cij <- matrix(0, nrow=nrow(cijr), ncol=ncol(cijr)); cij[upper.tri(cij)] <- cijr[upper.tri(cijr)];
rownames(cij) <- resnano_ras[resno_ras %in% inds_sca]
colnames(cij) <- resnano_ras[resno_ras %in% inds_sca]

# full ranked sca pairs
inds <- t(apply(as.matrix(sort(cij, decreasing=T)[1:(l*(l-1)/2)]), 1, function(x) { which(cij == x, arr.ind=T) }))
sca_ras <- as.data.frame(t(apply(inds, 1, function(x) { c(rownames(cij)[x],round(cij[x[1],x[2]],2)) })))
# add community
community_i <- names(membership_ras)[as.numeric(gsub("^.", "",sca_ras[,1]))]
community_j <- names(membership_ras)[as.numeric(gsub("^.", "",sca_ras[,2]))]
sca_ras <- cbind(sca_ras[,1:2], community_i, community_j, sca_ras[,3:ncol(sca_ras)])
colnames(sca_ras) <- c("a","b","community_i","community_j","sca")

# summary of inter and intra community pairs in top 100 pairs
tmp <- sca_ras[1:100, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(sca_ras[1:100, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_sca_ras_100 <- matrix(0, length(community_ras), length(community_ras))
rownames(tbl_sca_ras_100) <- community_ras; colnames(tbl_sca_ras_100) <- community_ras; 
tbl_sca_ras_100[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_sca_ras_100) <- diag(tbl_sca_ras_100)/2
tbl_sca_ras_100
         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8
a1/b1-b3        5  3  1  14     5  4  0  0  0
PL              3  3  0   0     0  2  0  0  0
SI              1  0  0   4     2  1  0  0  0
SII            14  0  4  22     5  9  0  0  0
b4-b6           5  0  2   5    11  2  0  1  3
a3              4  2  1   9     2  1  0  0  1
a4              0  0  0   0     0  0  0  0  0
a5              0  0  0   0     1  0  0  0  0
L8              0  0  0   0     3  1  0  0  1
sum(diag(tbl_sca_ras_100))
#[1] 43
sum(tbl_sca_ras_100[upper.tri(tbl_sca_ras_100)])
#[1] 57

# top 100 of eddm
tmp <- subseteddm_ras_gxp_beta0.5[1:100, c("community_j", "community_i")]
colnames(tmp) <- c("community_i", "community_j")
tmp <- rbind(subseteddm_ras_gxp_beta0.5[1:100, c("community_i", "community_j")], tmp)
tmp <- table(tmp)
tbl_eddm_ras_gxp_100 <- matrix(0, length(community_ras), length(community_ras))
rownames(tbl_eddm_ras_gxp_100) <- community_ras; colnames(tbl_eddm_ras_gxp_100) <- community_ras;
tbl_eddm_ras_gxp_100[rownames(tmp), colnames(tmp)] <- tmp
diag(tbl_eddm_ras_gxp_100) <- diag(tbl_eddm_ras_gxp_100)/2
tbl_eddm_ras_gxp_100
         a1/b1-b3 PL SI SII b4-b6 a3 a4 a5 L8
a1/b1-b3        5  1  6   7     0  0  0  1  0
PL              1  0  6  12     0  0  0  0  0
SI              6  6  5  23     0  0  0  0  0
SII             7 12 23  18     0 10  0  0  0
b4-b6           0  0  0   0     1  1  0  0  2
a3              0  0  0  10     1  1  1  0  0
a4              0  0  0   0     0  1  0  0  0
a5              1  0  0   0     0  0  0  0  0
L8              0  0  0   0     2  0  0  0  0
sum(diag(tbl_sca_ras_100))
#[1] 43
sum(tbl_sca_ras_100[upper.tri(tbl_sca_ras_100)])
#[1] 57

# top 100 pairs
pairs_eddm <- apply(as.matrix(subseteddm_ras_gxp_beta0.5[1:100,c("a","b")]),1,
  function(x) { paste0(x[1],"-",x[2]) }); names(pairs_eddm) <- 1:100
pairs_sca <- apply(as.matrix(sca_ras[1:100,c("a","b")]),1,
  function(x) { paste0(x[1],"-",x[2]) }); 

# common pairs found by sca and eddm
sum(pairs_sca %in% pairs_eddm)
#[1] 10
sca_ras[(1:100)[pairs_sca %in% pairs_eddm],]
subseteddm_ras_gxp_beta0.5[(1:100)[pairs_eddm %in% pairs_sca],-c(3:4)]

# smaller data.frame summarizing common pairs
inds <- match(pairs_eddm, pairs_sca)[!is.na(match(pairs_eddm, pairs_sca))]
sca_ras[(1:100)[inds],]
common_pairs <- cbind(subseteddm_ras_gxp_beta0.5[(1:100)[pairs_eddm %in% pairs_sca],c(1:2,5:9)],
  sca=sca_ras[(1:100)[inds],"sca"])

# fisher.test
fisher.test(matrix(c(10,90,90,l*(l-1)/2-100), 2,2))
#p-value = 6.517e-09

save(cijr, inds_sca, pairs_eddm, pairs_sca, l,
     sca_ras, tbl_sca_ras_100, tbl_eddm_ras_gxp_100, common_pairs, 
     file="sca_pairs.RData")



## name: res_eftu.r
## date: 11/18/2016

## Here I collect the residue name, number, membership of ras

load("/Users/hyangl/project/ras/results/2016/0715_ras_eftu/alignment/ali_raseftu.RData")

membership_eftu <- as.numeric(ali_raseftu["membership_eftu",ali_raseftu["membership_eftu",]!="0"])
community_eftu <- c("a1/b1-b3", "PL", "SI", "SII", "b4-b6", "a3", "a4", "a5", "L8","D2","D3")
for (i in 1:length(community_eftu)) { names(membership_eftu)[membership_eftu==i] <- community_eftu[i]}

# numeric vector and aa as name (e.g. "M" 72)
resno_eftu <- as.numeric(ali_raseftu["resno_eftu",ali_raseftu["membership_eftu",]!="0"])
names(resno_eftu) <- ali_raseftu["1TTT_A",ali_raseftu["membership_eftu",]!="0"]
# residue name + number format (e.g. "M72")
resnano_eftu <- paste0(names(resno_eftu), resno_eftu)

save(resno_eftu, resnano_eftu,
     membership_eftu, community_eftu,
     file = "res_eftu.RData")







## name: res_ras.r
## date: 11/20/2016

## Here I collect the residue name, number, membership of ras

load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")

membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
community_ras <- c("a1/b1-b3", "PL", "SI", "SII", "b4-b6", "a3", "a4", "a5", "L8")
for (i in 1:length(community_ras)) { names(membership_ras)[membership_ras==i] <- community_ras[i]}

# numeric vector and aa as name (e.g. "M" 72)
resno_ras <- as.numeric(ali["resno_ras",ali["membership_ras",]!="0"])
names(resno_ras) <- ali["1QRA_A",ali["membership_ras",]!="0"]
# residue name + number format (e.g. "M72")
resnano_ras <- paste0(names(resno_ras), resno_ras)

save(resno_ras, resnano_ras,
     membership_ras, community_ras,
     file = "res_ras.RData")







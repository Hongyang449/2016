## name: ali_raseftu.r
## date: 07/20/2016

load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
ali_rasgt <- ali 

seq_ras <- get.seq("5P21_A")
seq_eftu <- get.seq("1TTT_A")

# 405(resno in eftu)+42(gaps in ras)=448
ali <- matrix("0",nrow=7,ncol=447)
rownames(ali) <- c("1TTT_A", "5P21_A", "resno_eftu", "resno_ras",
  "membership_124","membership_eftu", "membership_ras")

# half-manually generate the alignment according to the PL_SI_SII fitting
ali[3,c(1:32,41:81,86:139,145:157,175:209,214:247,252:447)] <- 1:405
ali[4,c(9:40,68:71,73:144,146:174,201:213,236:251)] <- 1:166

# add residue names
ali[1,ali[3,]!=0] <- as.vector(seq_eftu$ali)
ali[2,ali[4,]!=0] <- as.vector(seq_ras$ali)
ali[1,ali[1,]==0] <- "-"
ali[2,ali[2,]==0] <- "-"

membership_ras <- ali_rasgt["membership_ras",ali_rasgt["membership_ras",]!=0]
ali["membership_ras",ali["resno_ras",]!=0] <- membership_ras
membership_eftu <- rep(0,405)
membership_eftu[c(9:17,26:38,66:80)] <- 1
membership_eftu[c(18:25)] <- 2
membership_eftu[c(39:65)] <- 3
membership_eftu[c(81:99)] <- 4
membership_eftu[c(100:109,130:136,169:175)] <- 5
membership_eftu[c(110:129)] <- 6
membership_eftu[c(145:168)] <- 7
membership_eftu[c(176:209)] <- 8
membership_eftu[c(137:144)] <- 9
membership_eftu[c(210:308)] <- 10
membership_eftu[c(309:405)] <- 11
ali["membership_eftu",which(ali["resno_eftu",]!="0")] <- membership_eftu

# membership_124
inds <- intersect(which(ali["membership_ras",]!=0),which(ali["membership_eftu",]!=0))
# !! 4 residues (117:120/137:140) of loop8 are different communities in two proteins
# ras 117:120 belong to b4-6 (I think this one is incorrect!!)
# eftu 137:140 belong to loop8
ali["membership_124", inds] <- ali["membership_ras",inds]

ali_raseftu <- ali

# write fasta file
write.fasta(alignment=list(ali=ali_raseftu[1:2,], id=rownames(ali_raseftu)[1:2]),
  ids=rownames(ali_raseftu)[1:2],file="fasta/5P21_A_1TTT_A.fa")

save(ali_raseftu, file="ali_raseftu.RData")



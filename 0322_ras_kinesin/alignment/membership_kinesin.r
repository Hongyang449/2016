## name: membership_kinesin.r
## date: 03/23/2016

## Here I try to partition kinesin into communities using ras network as a model
## The residue numbers of kinesin are from 4HNA_K.pdb

pdb_ki <- read.pdb("/Users/hyangl/project/ras/results/2016/0322_ras_kinesin/alignment/pdb/align_PL_SI_SII/aligned_4HNA_K.pdb")
pdb_ras <- read.pdb("/Users/hyangl/project/ras/results/2016/0322_ras_kinesin/alignment/pdb/align_PL_SI_SII/aligned_1QRA_A.pdb")
pdb_gt <- read.pdb("/Users/hyangl/project/ras/results/2016/0322_ras_kinesin/alignment/pdb/align_PL_SI_SII/aligned_1TAD_C.pdb")

# Find the seq of 4HNA; the last 5 inds are for ADP,MG,ALF, H2O*2
inds <- bounds(pdb_ki$atom[,"resno"], dup.inds=T)[,"start"][1:333]
seq_aa <- aa321(pdb_ki$atom[inds,"resid"])

# 1st row is residue number, 2nd row is membership, colname is aa type
membership_kinesin <- matrix(0, nrow=2, ncol=337)
colnames(membership_kinesin)[5:337] <- seq_aa
membership_kinesin[1,] <- 1:337
membership_kinesin[2,c(76:84,93:106,124:230)] <- 1
membership_kinesin[2,85:92] <- 2
membership_kinesin[2,190:204] <- 3
membership_kinesin[2,231:293] <- 4
membership_kinesin[2,c(294:304,5:17,49:56)] <- 5
membership_kinesin[2,305:325] <- 6
membership_kinesin[2,31:48] <- 7
membership_kinesin[2,57:75] <- 8
membership_kinesin[2,18:30] <- 9
membership_kinesin[2,326:337] <- 10
membership_kinesin[2,107:123] <- 11

membership_ki <- membership_kinesin[2,5:337]

save(membership_kinesin, file="membership_kinesin.RData")

# write a pdb in which the chain numbers are communities; the last 5 ADP,MG,ALF,H2O are chain A for now.
ch_aa <- vec2resno(c(LETTERS[membership_ki],rep("A",5)),pdb_ki$atom[,"resno"])
write.pdb(pdb_ki, chain=ch_aa, file="pdb/network/4HNA_network.pdb")

# write pdb for ras and gt
load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/comparison/ali.RData")
membership_ras <- as.numeric(ali["membership_ras",ali["membership_ras",]!="0"])
membership_gt <- as.numeric(ali["membership_gt",ali["membership_gt",]!="0"])

ch_aa <- vec2resno(c(LETTERS[membership_ras],rep("A",2)),pdb_ras$atom[,"resno"])
write.pdb(pdb_ras, chain=ch_aa, file="pdb/network/1QRA_network.pdb")

attach(transducin)
inds <- pdbs$resno[1,!is.gap(pdbs)]
pdb_305 <- trim.pdb(pdb_gt, resno=inds)
ch_aa <- vec2resno(c(LETTERS[membership_gt]),pdb_305$atom[,"resno"])
write.pdb(pdb_305, chain=ch_aa, file="pdb/network/1TAD_network.pdb")


# community 1 (b1-3,a1)
# 1:9, 18:24, 38:56
# 76:84, 93:106, 124:230

# community 2 (PL)
# 10:17
# 85:92

# community 3 (SI)
# 25:37
# 190:204

# community 4 (SII)
# 57:75
# 231:293

# community 5 (b4-6)
# 76:85, 110-121, 140:148
# 294:304, 5:17, 49:56

# community 6 (SIII/a3)
# 86:109
# 305:325

# community 7 (a4) (3beta in kinesin)
# 131:139
# 31:48

# community 8 (a5)
# 149:166
# 57:75

# community 9 (loop8)
# 122:130
# 18:30

# community 10 (NL)
# NA
# 326:337

# community 11 (loop3?) (a2 in kinesin)
# NA
# 107:123





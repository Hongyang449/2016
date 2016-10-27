## name: 1101_ras_kinesin 
## date: 10/26/2015

# prepare ras/kinesin pdb
pdb <- get.pdb("4HNA", path="pdb/raw")
pdbsplit(pdb, path="pdb/split_chain")

cp /Users/hyangl/project/ras/results/2015/0730_ras_transducin/alignment/aligned_1QRA_A.pdb pdb/
cp pdb/split_chain/4HNA_K.pdb pdb/4HNA_K.pdb

# structural alignment (generate results.html results.pdb)
mustang -i pdb/aligned_1QRA_A.pdb pdb/4HNA_K.pdb -F 'fasta' -r 'ON' -o '1QRA_A_4HNA_K'

# R trim.pdb
x <- read.pdb("1QRA_A_4HNA_K.pdb")
write.pdb(trim.pdb(x, chain="A"), file="aligned_1QRA_A.pdb")
write.pdb(trim.pdb(x, chain="B"), file="aligned_4HNA_K.pdb")

## The structural alignment of ras and kinesin is terrible!
## Let us try structural superposition - superpose beta1, p-loop
pdb_ras <- read.pdb("pdb/align_PL_B1/aligned_1QRA_A.pdb")
pdb_gt <- read.pdb("pdb/align_PL_B1/aligned_1TAD_C.pdb")
pdb_ki <- read.pdb("pdb/split_chain/4HNA_K.pdb")

# Here I only superpose beta1 and P-loop
inds_ras <- atom.select(pdb_ras, resno=1:17, elety="CA")
inds_ki <- atom.select(pdb_ki, resno=76:92, elety="CA")

pdb_ki_ali <- pdb_ki
pdb_ki_ali$xyz <- fit.xyz(pdb_ras$xyz, pdb_ki$xyz, inds_ras$xyz, inds_ki$xyz)
pdb_ki_ali$atom[,c("x","y","z")] <- t(matrix(pdb_ki_ali$xyz[1,], nrow=3))

write.pdb(pdb_ki_ali, file="pdb/align_PL_B1/aligned_4HNA_K.pdb")


## date: 03/23/2016
## Here I use conserved residues in PL,SI,SII to align ras and ki
## Actually the alignment above using PL,B1 is pretty good!

## Five G domain motifs
# G1. P-loop(Walker A): G-x(4)-GK-[TS] 10:17
# G2. Switch I: T 35
# G3. Switch II: D-xx-G 57:60
# G4. specific guanine binding: NK-x-DL 116:120
# G5. specific guanine binding: SAK 145:147

pdb_ras <- read.pdb("pdb/align_PL_SI_SII/aligned_1QRA_A.pdb")
pdb_gt <- read.pdb("pdb/align_PL_SI_SII/aligned_1TAD_C.pdb")
pdb_ki <- read.pdb("pdb/split_chain/4HNA_K.pdb")

# conserved GTP/ATP binding sites in PL,SI,SII
inds_ras <- atom.select(pdb_ras, resno=c(10:17,35,57:60), elety="CA")
inds_ki <- atom.select(pdb_ki, resno=c(85:92,202,231:234), elety="CA")

pdb_ki_ali <- pdb_ki
pdb_ki_ali$xyz <- fit.xyz(pdb_ras$xyz, pdb_ki$xyz, inds_ras$xyz, inds_ki$xyz)
pdb_ki_ali$atom[,c("x","y","z")] <- t(matrix(pdb_ki_ali$xyz[1,], nrow=3))

write.pdb(pdb_ki_ali, file="pdb/align_PL_SI_SII/aligned_4HNA_K.pdb")















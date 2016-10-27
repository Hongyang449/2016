## name: 0715_ras_eftu
## date: 07/15/2016

# prepare ras/eftu pdb (GTP bound)
pdb <- get.pdb("1TTT", path="pdb/raw")
pdbsplit(pdb, path="pdb/split_chain")
pdb <- get.pdb("1TUI", path="pdb/raw")
pdbsplit(pdb, path="pdb/split_chain")

cp /Users/hyangl/project/ras/results/2015/0730_ras_transducin/alignment/aligned_1QRA_A.pdb pdb/
cp pdb/split_chain/1TTT_A.pdb pdb/1TTT_A.pdb

# structural alignment (generate results.html results.pdb)
# cd mustang_ali
mustang -i ../pdb/aligned_1QRA_A.pdb ../pdb/1TTT_A.pdb -F 'fasta' -r 'ON' -o '1QRA_A_1TTT_A'

# R trim.pdb
x <- read.pdb("1QRA_A_1TTT_A.pdb")
write.pdb(trim.pdb(x, chain="A"), file="aligned_1QRA_A.pdb")
write.pdb(trim.pdb(x, chain="B"), file="aligned_1TTT_A.pdb")

## date: 07/18/2016

## The structural alignment of ras and eftu is ok, but..
## Here I use conserved residues in PL,SI,SII to align ras and eftu

## Five G domain motifs
# G1. P-loop(Walker A): G-x(4)-GK-[TS] 10:17
# G2. Switch I: T 35
# G3. Switch II: D-xx-G 57:60
# G4. specific guanine binding: NK-x-DL 116:120
# G5. specific guanine binding: SAK 145:147

## eftu
# G1. G-x(4)-GKT 18:25
# G2. T 62
# G3. D-xx-G 81:84 
# G4. NK-x-DM 136:140
# G5. SAL 174:176

# cp /Users/hyangl/project/ras/results/2016/0503_ras_transducin/alignment/mustang_ali/aligned_5P21_A.pdb pdb/align_PL_SI_SII/
pdb_ras <- read.pdb("pdb/align_PL_SI_SII/aligned_5P21_A.pdb")
pdb_eftu <- read.pdb("pdb/split_chain/1TTT_A.pdb")

# conserved GTP/ATP binding sites in PL,SI,SII
inds_ras <- atom.select(pdb_ras, resno=c(10:17,35,57:60), elety="CA")
inds_eftu <- atom.select(pdb_eftu, resno=c(18:25,62,81:84), elety="CA")

pdb_eftu_ali <- pdb_eftu
pdb_eftu_ali$xyz <- fit.xyz(pdb_ras$xyz, pdb_eftu$xyz, inds_ras$xyz, inds_eftu$xyz)
pdb_eftu_ali$atom[,c("x","y","z")] <- t(matrix(pdb_eftu_ali$xyz[1,], nrow=3))

write.pdb(pdb_eftu_ali, file="pdb/align_PL_SI_SII/aligned_1TTT_A.pdb")

# fit 1TUI_A to 1TTT BUT excluding SI
pdb_gtp <- read.pdb("pdb/align_PL_SI_SII/aligned_1TTT_A.pdb")
pdb_gdp <- read.pdb("pdb/split_chain/1TUI_A.pdb")

# conserved GTP/ATP binding sites in PL,SI,SII
inds_gtp <- atom.select(pdb_gtp, resno=c(9:38,66:209), elety="CA")
inds_gdp <- atom.select(pdb_gdp, resno=c(9:38,66:209), elety="CA")

pdb_ali <- pdb_gdp
pdb_ali$xyz <- fit.xyz(pdb_gtp$xyz, pdb_gdp$xyz, inds_gtp$xyz, inds_gdp$xyz)
pdb_ali$atom[,c("x","y","z")] <- t(matrix(pdb_ali$xyz[1,], nrow=3))

write.pdb(pdb_ali, file="pdb/align_PL_SI_SII/aligned_1TUI_A.pdb")


## date: 07/19/2016

# rum vmd drawing lines
vmd -m pdb/align_PL_SI_SII/*.pdb -e draw_lines.tcl



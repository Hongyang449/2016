## name: draw_lines_between_atoms.r
## date: 07/19/2016

# Here I want to draw lines between closest CAs according to align_PL_SI_SII 
# Thus I can align them residue-wise and find residue/community counterparts
# vmd -m pdb/align_PL_SI_SII/*.pdb -e draw_lines.tcl


pdb_ras <- read.pdb("/Users/hyangl/project/ras/results/2016/0715_ras_eftu/alignment/pdb/align_PL_SI_SII/aligned_1QRA_A.pdb")
pdb_eftu <- read.pdb("/Users/hyangl/project/ras/results/2016/0715_ras_eftu/alignment/pdb/align_PL_SI_SII/aligned_1QRA_A.pdb")

# xyz of ca
xyz_ras <- pdb_ras$atom[pdb_ras$calpha,c("x","y","z")]
xyz_eftu <- pdb_eftu$atom[pdb_eftu$calpha,c("x","y","z")]

# ca distance matrix between ras and eftu (166 * 405)
mat <- dist.xyz(as.vector(t(xyz_ras)),as.vector(t(xyz_eftu)))

# inds is the closest residues of eftu to 1:166 residue of ras
inds <- apply(mat,1,order)[1,]

# cat the script!
out <- "draw_lines.tcl"
cat("", file=out)
for (i in 1:length(inds)) {
  cat("graphics top line { ", xyz_ras[i,1], " ", xyz_ras[i,2], " ", xyz_ras[i,3],
  " } { ", xyz_eftu[inds[i],1], " ", xyz_eftu[inds[i],2], " ", xyz_eftu[inds[i],3],
  " } width 5 style solid\n", sep="", file=out, append=TRUE)
}





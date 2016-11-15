## ./pdb/
# create links to aligned pdb
ln -s /Users/hyangl/project/ras/results/2016/0503_ras_transducin/alignment/mustang_ali/aligned_5P21_A.pdb aligned_5P21_A.pdb
ln -s /Users/hyangl/project/ras/results/2016/0503_ras_transducin/alignment/mustang_ali/aligned_1TND_B.pdb aligned_1TND_B.pdb
ln -s /Users/hyangl/project/ras/results/2016/0715_ras_eftu/alignment/pdb/align_PL_SI_SII/aligned_1TTT_A.pdb aligned_1TTT_A.pdb
ln -s /Users/hyangl/project/ras/results/2016/0715_ras_eftu/alignment/pdb/align_PL_SI_SII/aligned_1TUI_A.pdb aligned_1TUI_A.pdb

# load ali and find aligned residues (eftu add a4)
# presented as "Cartoon"
paste(ali["resno_ras",ali["membership_146",]!="0"], collapse=" ")
paste(ali["resno_gt",ali["membership_146",]!="0"], collapse=" ")
paste(c(ali_raseftu["resno_eftu",ali_raseftu["membership_124",]!="0"], 144:161), collapse=" ")
# presented as "Tube"
paste(ali["resno_ras",(ali["membership_146",]=="0"&ali["resno_ras",]!="0")], collapse=" ")
paste(ali["resno_gt",(ali["membership_146",]=="0"&ali["resno_gt",]!="0")], collapse=" ")
inds <- which(!(ali_raseftu["resno_eftu",(ali_raseftu["membership_124",]=="0"&
  ali_raseftu["resno_eftu",]!="0")] %in% 144:161))
paste(ali_raseftu["resno_eftu",(ali_raseftu["membership_124",]=="0"&ali_raseftu["resno_eftu",]!="0")][inds],
  collapse=" ")
# save vmd session as pdb/aligned_ras_gt_eftu.vmd





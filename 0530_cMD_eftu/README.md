## name: cMD_mutations
## date: 05/30/2016

# wt
# GTP(1TTT_A) 1:405
# GDP(1TUI_A) 9:405 !!
cp /Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pdbs/corefit/1TTT_A.pdb_flsq.pdb pdb/wt/1TTT.pdb
cp /Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/pdbs/corefit/1TUI_A.pdb_flsq.pdb pdb/wt/1TUI.pdb

# generate pdb/mutation/1TTT_I93A.pdb using pymol
pymol pdb/split_chain/1TTT_A.pdb
select I93, resi 93
# Wizard -> mutagenesis
# pick a residue (I93) -> select mutation to Ala ->
# select conformation with no/small clash -> apply -> done ->
# save molecule pdb/mutation/1TTT_I93.pdb


cp 2015/1127_cMD_Ras/sample/5P21_R164A_setup_cMD.sh sample/


cp /Users/hyangl/project/ras/results/2015/0206_cMD/dccm/*.r dccm/



sed -e 's/R75A/V128A/g' < sample/1TTT_R75A_setup_cMD.sh > sample/1TTT_V128A_setup_cMD.sh

sed -e 's/D207/R75/g' < dccm_1TTT_D207A/pdb/README.md > dccm_1TTT_R75A/pdb/README.md
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/write_ncdf.r > dccm_1TTT_R75A/write_ncdf.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/cmap.r > dccm_1TTT_R75A/cmap.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/cij_rasd.r > dccm_1TTT_R75A/cij_rasd.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/rmsd_rmsf.r > dccm_1TTT_R75A/rmsd_rmsf.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/rmsd_plot.r > dccm_1TTT_R75A/rmsd_plot.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/rmsf_plot_mean.r > dccm_1TTT_R75A/rmsf_plot_mean.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/dccm_plot.r > dccm_1TTT_R75A/dccm_plot.r
sed -e 's/D207/R75/g' < dccm_1TTT_D207A/net_eftu_rasd_2f.r > dccm_1TTT_R75A/net_eftu_2f_wilcox.r





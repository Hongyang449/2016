## name: write_pdbs_vmd.r
## date: 05/30/2016

## Here I want to generate a .vmd file to present colored superposed pdbs
## vmd -e XXX.vmd

####################
## colored by PCA ##
####################

out <- "vmd/eftu_pdbs_2_grps.vmd"
pwd <- "/Users/hyangl/project/ras/results/2016/0415_pdbs/eftu/"

load(paste0(pwd, "pca.RData"))

vcolors <- rep(0, length(pdbs_eftu_new$id))
#vcolors[as.character(eftu_annotate[,"colorLig"]=="blue"] <- 11
#vcolors[as.character(eftu_annotate[,"colorLig"]=="orange"] <- 3
vcolors[as.character(eftu_annotate[,"colorLig"])=="red"] <- 1
vcolors[as.character(eftu_annotate[,"colorLig"])=="green"] <- 7

## write .vmd
cat("", file=out)
for(i in 1:length(pdbs_eftu_new$id)) {
  filename <- paste0(pwd, "pdbs/corefit/", pdbs_eftu_new$id[i], ".pdb_flsq.pdb")
  resno <- pdbs_eftu_new$resno[i, gaps.res$f.inds]
  cat("mol new ", filename, "\n",
    "mol delrep 0 top\n",
    "mol selection \"resid ", resno[1], " to ", resno[length(resno)],"\"\n",
    "mol color ColorID ", vcolors[i], "\n",
    "mol rep newcartoon\n",
    "mol addrep top\n", sep="", file=out, append=TRUE)
}
cat("color Display {Background} white\n",
    "display depthcue off\n", sep="", file=out, append=TRUE)

#########################
## colored identically ##
#########################

out <- "vmd/eftu_pdbs.vmd"

cat("", file=out)
for(i in 1:length(pdbs_eftu_new$id)) {
  filename <- paste0(pwd, "pdbs/corefit/", pdbs_eftu_new$id[i], ".pdb_flsq.pdb")
  resno <- pdbs_eftu_new$resno[i, gaps.res$f.inds]
  cat("mol new ", filename, "\n",
    "mol delrep 0 top\n",
    "mol selection \"resid ", resno[1], " to ", resno[length(resno)],"\"\n",
    "mol color ColorID 0\n",
    "mol rep newcartoon\n",
    "mol addrep top\n", sep="", file=out, append=TRUE)
}
cat("color Display {Background} white\n",
    "display depthcue off\n", sep="", file=out, append=TRUE)













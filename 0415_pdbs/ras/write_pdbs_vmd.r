## name: write_pdbs_vmd.r
## date: 04/15/2016

## Here I want to generate a .vmd file to present colored superposed pdbs
## vmd -e XXX.vmd

####################
## colored by PCA ##
####################

out <- "vmd/ras_pdbs_4_grps.vmd"
pwd <- "/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/"

load(paste0(pwd, "pdbs.RData"))

vcolors <- rep(0, length(pdbs_ras_new$id))
vcolors[cols=="purple"] <- 11
vcolors[cols=="red"] <- 1
vcolors[cols=="orange"] <- 3
vcolors[cols=="green"] <- 7

gaps.res <- gap.inspect(pdbs_ras_new$ali)

## write .vmd
cat("", file=out)
for(i in 1:length(pdbs_ras_new$id)) {
  filename <- paste0(pwd, "corefit/", pdbs_ras_new$id[i], ".pdb_flsq.pdb")
  resno <- pdbs_ras_new$resno[i, gaps.res$f.inds]
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

out <- "vmd/ras_pdbs.vmd"

cat("", file=out)
for(i in 1:length(pdbs_ras_new$id)) {
  filename <- paste0(pwd, "corefit/", pdbs_ras_new$id[i], ".pdb_flsq.pdb")
  resno <- pdbs_ras_new$resno[i, gaps.res$f.inds]
  cat("mol new ", filename, "\n",
    "mol delrep 0 top\n",
    "mol selection \"resid ", resno[1], " to ", resno[length(resno)],"\"\n",
    "mol color ColorID 10\n",
    "mol rep newcartoon\n",
    "mol addrep top\n", sep="", file=out, append=TRUE)
}
cat("color Display {Background} white\n",
    "display depthcue off\n", sep="", file=out, append=TRUE)













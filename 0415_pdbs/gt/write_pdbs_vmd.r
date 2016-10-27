## name: write_pdbs_vmd.r
## date: 04/15/2016

## Here I want to generate a .vmd file to present colored superposed pdbs
## vmd -e XXX.vmd

####################
## colored by PCA ##
####################

out <- "vmd/gt_pdbs_4_grps.vmd"
pwd <- "/Users/hyangl/project/ras/results/2016/0415_pdbs/gt/"

attach(transducin)

vcolors <- rep(0, length(pdbs$id))
vcolors[annotation[,"state3"]=="GDI"] <- 0
vcolors[annotation[,"state3"]=="GTP"] <- 1
vcolors[annotation[,"state3"]=="GDP"] <- 7

gaps.res <- gap.inspect(pdbs$ali)

## write .vmd
cat("", file=out)
for(i in 1:length(pdbs$id)) {
  filename <- paste0(pwd, "corefit/", pdbs$id[i], ".pdb_flsq.pdb")
  resno <- pdbs$resno[i, gaps.res$f.inds]
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

out <- "vmd/gt_pdbs.vmd"

cat("", file=out)
for(i in 1:length(pdbs$id)) {
  filename <- paste0(pwd, "corefit/", pdbs$id[i], ".pdb_flsq.pdb")
  resno <- pdbs$resno[i, gaps.res$f.inds]
  cat("mol new ", filename, "\n",
    "mol delrep 0 top\n",
    "mol selection \"resid ", resno[1], " to ", resno[length(resno)],"\"\n",
    "mol color ColorID 10\n",
    "mol rep newcartoon\n",
    "mol addrep top\n", sep="", file=out, append=TRUE)
}
cat("color Display {Background} white\n",
    "display depthcue off\n", sep="", file=out, append=TRUE)













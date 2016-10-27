## name: write_pdbs.r
## date: 04/16/2016

## Here I write fitted .pdb based on coordinates of attached pdbs

attach(transducin)
## gt 53 chains from 36 structures

files <- get.pdb(pdbs$id, path="pdbs/")
pdbsplit(files, path="pdbs/split_chain/")
# add path for pdbs$id then you can pdbfit() to write fitted .pdb
pdbs$id <- paste("pdbs/split_chain/", pdbs$id, ".pdb", sep="")
pdbfit(pdbs, outpath="corefit")

## You can corefit again
#files <- paste("pdbs/split_chain/", pdbs$id, ".pdb", sep="")
#pdbs1 <- pdbaln(files, ncore=10, outfile="fasta/gt_pdbs.fa")
#core <- core.find(pdbs1)
#core.inds <- print(core)
#xyz <- pdbfit(pdbs1, inds=core.inds, outpath="corefit/")




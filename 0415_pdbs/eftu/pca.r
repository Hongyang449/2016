## name: pca.r
## date: 05/27/2016

source("/Users/hyangl/project/ras/results/2015/functions/ligand_functions.R")
anno.terms <- c("structureId","chainId","compound","ligandId","source",
  "resolution","experimentalTechnique","publicationYear","title")

# Extract sequence and BLAST against PDB database
seq_1d8t <- pdbseq( read.pdb("1D8T") )
blast <- blast.pdb(seq_1d8t)
# [1] 608 hits (some redundance)
length(unique(substr(blast$pdb.id,1,4)))
# [1] 212 unique PDB structures

# Select subset of hits based on the Evalues of BLAST (here it is cut between e-103 and e-66)
plot(blast)
#  * Possible cutoff values:    709 15
#            Yielding Nhits:    224 608
#
#  * Chosen cutoff value of:    709
#            Yielding Nhits:    224
hits <- plot(blast, cutoff=709)
unq.ids <- unique(substr(hits$pdb.id,1,4))
length(hits$pdb.id)
# [1] 224
length(unq.ids)
# [1] 77

# Download these PDB files and split them into seperate chains
files <- get.pdb(hits$pdb.id, path="pdbs/")
pdbsplit(files, path="pdbs/split_chain/")

# Align the split chains
files <- paste("pdbs/split_chain/", hits$pdb.id, ".pdb", sep="")
pdbs_eftu <- pdbaln(files, ncore=10, outfile="fasta/eftu_raw.fa")
eftu_annotate_raw <- pdb.annotate(unique(substr(basename(pdbs_eftu$id),1,4)), anno.terms=anno.terms)
write.csv(eftu_annotate_raw,file="annotation/eftu_raw.csv")
# cannot retrieve some PDBs:
y <- unique(substr(basename(eftu_annotate_raw[,1]),1,4))
x <- unique(substr(basename(pdbs_eftu$id),1,4))
x[which(is.na(match(x,y)))]
# [1] "3FIH" "3IZV" "3IZW" "4ABR" "2Y0U" "2Y0W" "3FIC" "2XQD" "2WRN" "2WRQ"
#[11] "2Y0Y" "2Y10" "2Y12" "2Y14" "2Y16" "2Y18"
# These PDBs are renamed, for now we skip them. If needed, come back later!!

# rm EM strucutres
inds <- which(eftu_annotate_raw[,"experimentalTechnique"] == "X-RAY DIFFRACTION")
# be careful to find ids; sort ids; unique to rm duplicates (a lot of duplicates even in hits$pdb.id)
ids <- unique(sort(hits$pdb.id[!is.na(match(substr(hits$pdb.id,1,4),
  unique(eftu_annotate_raw[inds,"structureId"])))]))
files <- paste("pdbs/split_chain/", ids, ".pdb", sep="")
pdbs_eftu <- pdbaln(files=files, outfile="fasta/eftu_xray0.fa", ncore=10)
eftu_annotate_xray0 <- pdb.annotate(unique(substr(basename(pdbs_eftu$id),1,4)), anno.terms=anno.terms)
write.csv(eftu_annotate_xray0,file="annotation/eftu_xray0.csv")

# seaview eftu_xray.fa
# rm sequences with long gaps (and missing PL motif)
gaps.res <- gap.inspect(pdbs_eftu$ali)
inds1 <- which(!gaps.res$row >= 911)
files <- pdbs_eftu$id[inds1]
pdbs_eftu <- pdbaln(files=files, outfile="fasta/eftu_xray1.fa", ncore=10)
eftu_annotate_xray1 <- pdb.annotate(unique(substr(basename(pdbs_eftu$id),1,4)), anno.terms=anno.terms)
write.csv(eftu_annotate_xray1,file="annotation/eftu_xray1.csv")

# rm PDBs with gaps 
ids_rm <- c("3U6B_A","3U6B_B","4IW3_B","4IW3_K")
inds2 <- !(substring(basename(pdbs_eftu$id),1,6) %in% ids_rm)
files <- pdbs_eftu$id[inds2]
pdbs_eftu <- pdbaln(files=files, outfile="fasta/eftu_xray2.fa", ncore=10)
eftu_annotate_xray2 <- pdb.annotate(unique(substr(basename(pdbs_eftu$id),1,4)), anno.terms=anno.terms)
write.csv(eftu_annotate_xray2,file="annotation/eftu_xray2.csv")

# seaview eftu_xray2.fa
# manually rm PDBs with missing residues (mainly in SI) result in 34 chains from 23 PDBs
pdbs_eftu <- read.fasta("fasta/eftu_xray_tmp.fa")
# pdbaln again
files <- sort(pdbs_eftu$id)
pdbs_eftu <- pdbaln(files=files, outfile="fasta/eftu.fa", ncore=10)
eftu_annotate <- pdb.annotate(unique(substr(basename(pdbs_eftu$id),1,4)), anno.terms=anno.terms)
# add lig annotation
inds3 <- which(rownames(eftu_annotate) %in% substring(basename(pdbs_eftu$id),1,6))
ligs_raw <- eftu_annotate[inds3, "ligandId"]
ligs <- unlist(lapply(as.list(1:length(ligs_raw)), function(x) {
  tmp <- class.convert.yao(each.lig(ligs_raw[x]))
  # if NA or APO, it returns character(0) and unlist() will skip that row 
  if(identical(tmp, character(0))) return("APO")
  else return(tmp)
}))
# add lig color
col_lig <- rep("gray", length(ligs))
col_lig[ligs=="TP"] <- "red"
col_lig[ligs=="DP"] <- "green"
col_lig[ligs=="APO"] <- "black"
# add source color
col_source <- rep("gray", length(ligs))
tmp <- table(eftu_annotate[inds3,"source"])
for (i in 1:length(ligs)) {
  col_source[i] <- vmd.colors()[which(eftu_annotate[inds3[i],"source"]==names(tmp))]
}
eftu_annotate <- cbind(eftu_annotate[inds3,], "ligs"=ligs, "colorLig"=col_lig, "colorSource"=col_source)
# this annotation table only contain eftu chains no other chains
write.csv(eftu_annotate,file="annotation/eftu.csv")

####################

# PCA finally
# find core
core <- core.find(pdbs_eftu)
# Select core residues based on volume-length plot; default <= 1 A^3
core.inds <- print(core)
# 91 positions
# Superpose the structures on the core residues
xyz <- pdbfit(pdbs_eftu, inds=core.inds, outpath="pdbs/corefit/")
# Write core structure
write.pdb(xyz = xyz[1, core.inds$xyz], resno=pdbs_eftu$resno[1, core.inds$atom], file = "eftu_core.pdb")
# assign fitted xyz to pdbs
pdbs_eftu$xyz <- xyz

pca_eftu <- pca(xyz, rm.gaps=TRUE)
# pca_eftu <- pca(pdbs_eftu$xyz, rm.gaps=TRUE)
# pca_eftu <- pca.xyz(xyz[, gaps.pos$f.inds])

plot.pca(pca_eftu, col=as.character(eftu_annotate[,"colorLig"]))
dev.copy2pdf(file="figures/pca_eftu_ligs.pdf")
plot.pca(pca_eftu, col=as.character(eftu_annotate[,"colorSource"]))
dev.copy2pdf(file="figures/pca_eftu_source.pdf")

# PCA on all structures
gaps.res <- gap.inspect(pdbs_eftu$ali)
gaps.pos <- gap.inspect(pdbs_eftu$xyz)

# Write PC trajectory
a <- mktrj.pca(pca_eftu, pc=1, file="pca/eftu_pc1.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )
b <- mktrj.pca(pca_eftu, pc=2, file="pca/eftu_pc2.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )
c <- mktrj.pca(pca_eftu, pc=3, file="pca/eftu_pc3.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )

# new pdbs
source("/Users/hyangl/project/ras/results/2015/functions/pdbaln_new.R")
pdbs_eftu_new <- pdbaln_new(files=pdbs_eftu$id, outfile=NULL, ncore=1)

eftu_annotate <- as.matrix(eftu_annotate)

save(blast, core, core.inds, gaps.pos, gaps.res, 
     pca_eftu, pdbs_eftu, pdbs_eftu_new, eftu_annotate,
     file="pca.RData")

save(blast, pdbs_eftu, pdbs_eftu_new,
     file="pdbs.RData")

save(eftu_annotate,
     file="eftu_annotate.RData")


##############################
##############################

## I also tried to use eftu_xray2 PDBs (76 chains but missing SI)
## mainly GTP/GDP/APO states can be distinguished but some outliers are mixed (1HA3,4PC2/3/6,4FWT,3VNU)
pdbs_eftu <- read.fasta("fasta/eftu_xray2.fa")
files <- sort(pdbs_eftu$id)
pdbs_eftu <- pdbaln(files=files, outfile=NULL, ncore=10)
eftu_annotate <- pdb.annotate(unique(substr(basename(pdbs_eftu$id),1,4)), anno.terms=anno.terms)
# add lig annotation
inds3 <- which(rownames(eftu_annotate) %in% substring(basename(pdbs_eftu$id),1,6))
ligs_raw <- eftu_annotate[inds3, "ligandId"]
ligs <- unlist(lapply(as.list(1:length(ligs_raw)), function(x) {
  tmp <- class.convert.yao(each.lig(ligs_raw[x]))
  # if NA or APO, it returns character(0) and unlist() will skip that row 
  if(identical(tmp, character(0))) return("APO")
  else return(tmp)
}))
# add lig color
col_lig <- rep("gray", length(ligs))
col_lig[ligs=="TP"] <- "red"
col_lig[ligs=="DP"] <- "green"
col_lig[ligs=="APO"] <- "black"
# add source color
col_source <- rep("gray", length(ligs))
tmp <- table(eftu_annotate[inds3,"source"])
for (i in 1:length(ligs)) {
  col_source[i] <- vmd.colors()[which(eftu_annotate[inds3[i],"source"]==names(tmp))]
}
eftu_annotate <- cbind(eftu_annotate[inds3,], "ligs"=ligs, "colorLig"=col_lig, "colorSource"=col_source)

# PCA finally
# find core
core <- core.find(pdbs_eftu)
# Select core residues based on volume-length plot; default <= 1 A^3
core.inds <- print(core)
# Superpose the structures on the core residues
xyz <- pdbfit(pdbs_eftu, inds=core.inds, outpath=NULL)

pca_eftu <- pca(xyz, rm.gaps=TRUE)
# pca_eftu <- pca(pdbs_eftu$xyz, rm.gaps=TRUE)
# pca_eftu <- pca.xyz(xyz[, gaps.pos$f.inds])

plot.pca(pca_eftu, col=as.character(eftu_annotate[,"colorLig"]))
dev.copy2pdf(file="figures/xray2_pca_eftu_ligs.pdf")
plot.pca(pca_eftu, col=as.character(eftu_annotate[,"colorSource"]))

# PCA on all structures
gaps.res <- gap.inspect(pdbs_eftu$ali)
gaps.pos <- gap.inspect(pdbs_eftu$xyz)

# Write PC trajectory
a <- mktrj.pca(pca_eftu, pc=1, file="eftu_pc1.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )
b <- mktrj.pca(pca_eftu, pc=2, file="eftu_pc2.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )
c <- mktrj.pca(pca_eftu             , pc=3, file="eftu_pc3.pdb",
               resno = pdbs_eftu$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_eftu$ali[1, gaps.res$f.inds]) )

plot(cbind(pca_eftu$z[,1],pca_eftu$z[,2]), cex=1.0, col=as.character(eftu_annotate[,"colorLig"]), 
  xlab="pc1", ylab="pc2", main="pca_eftu_pc12_ligs")
identify(pca_eftu$z[,1],pca_eftu$z[,2],labels=substr(basename(pdbs_eftu$id),1,6))
dev.copy2pdf(file="figures/pc12_eftu_ligs.pdf")







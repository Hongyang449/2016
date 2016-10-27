## name: ras_pdbs.r
## date: 09/06/2016

##########
## pdbs ##
##########

# extract sequence and BLAST against PDB database
seq_121p <- pdbseq( read.pdb("121p") )
blast <- blast.pdb(seq_121p)
# [1] 319 hits
length(unique(substr(blast$pdb.id,1,4)))
# [1] 240 unique PDB structures

# Select subset of hits based on the Evalues of BLAST
hits <- plot(blast, cutoff=233)
ids <- unique(substr(hits$pdb.id,1,4))
length(hits$pdb.id)
# [1] 233
length(ids)
# [1] 179

# fetch and align raw hits
files <- get.pdb(hits$pdb.id, path="pdbs/")
pdbsplit(files, path="pdbs/split_chain/")
files <- paste("pdbs/split_chain/", hits$pdb.id, ".pdb", sep="")
pdbs <- pdbaln(files, ncore=8, outfile="fasta/ras_raw.fa")

# index for 1:166
index1 <- pdbs$ali["pdbs/split_chain/5P21_A.pdb",] != "-"

# rm pdbs with gaps
gaps <- gap.inspect(pdbs$ali)
index2 <- apply(!gaps$bin[,index1], 1, all)
files <- files[index2]
pdbs <- pdbaln(files=files, ncore=8, outfile=NULL)

# pdb annotation
ras_annotate_all <- pdb.annotate(unique(substr(basename(pdbs$id),1,4)), anno.terms=
      c("structureId", "chainId", "resolution", "ligandId", "ligandName", "compound", "source",
        "scopDomain", "spaceGroup", "title", "experimentalTechnique"))
ras_annotate_all <- ras_annotate[order(rownames(ras_annotate)),]

# remove NMR structures
rmindex <- NULL
ids_nmr <- rownames(ras_annotate)[(ras_annotate[,"experimentalTechnique"]=="SOLUTION NMR")]
for (i in ids_nmr) {
  rmindex <- c(grep(i, pdbs$id), rmindex)
}
files <- sort(pdbs$id[-rmindex])
pdbs <- pdbaln(files=files, outfile=NULL, ncore=10)
length(pdbs$id)
# [1] 143
length(unique(substr(basename(pdbs$id),1,4)))
# [1] 121

# pdbs_new rename pdbs as 5P21_A
source("/Users/hyangl/project/ras/results/2015/functions/pdbaln_new.R")
pdbs_new <- pdbaln_new(files=files, outfile="fasta/ras_xray.fa", ncore=1)

pdbs_ras <- pdbs
pdbs_ras_new <- pdbs_new
save(blast, pdbs_ras, pdbs_ras_new, cols,
     file="pdbs.RData")

####################
## pdb annotation ##
####################

ras_annotate <- pdb.annotate(unique(substr(basename(files),1,4)), anno.terms=
      c("structureId", "chainId", "resolution", "ligandId", "ligandName", "compound", "source",
        "scopDomain", "spaceGroup", "title", "experimentalTechnique"))
ras_annotate <- ras_annotate[order(rownames(ras_annotate)),]
ras_annotate_trim <- ras_annotate[rownames(ras_annotate) %in% pdbs_ras_new$id,]
save(ras_annotate, ras_annotate_all, ras_annotate_trim,
  file="tmp_ras_annotate.RData")

source("/Users/hyangl/project/ras/results/2015/functions/ligand_functions.R")

cl <- rep(NA, nrow(ras_annotate_trim))
for(i in 1:nrow(ras_annotate_trim)) {
  tmp <- class.convert.yao( each.lig(ras_annotate_trim[i,"ligandId"]) )
  if(length(tmp) > 0) {
    if (length(tmp) > 1) {
      print( paste(" Row:",i))
      print(tmp)
      tmp <- paste(unique(tmp), collapse=" ")
    }
    cl[i] <- tmp
  }
}
# transition state
cl[which(nchar(cl) > 2)] = "TS"
names(cl) = rownames(ras_annotate_trim)

cols <- cl
cols[cl=="TP"] <- "red"
cols[cl=="DP"] <- "green"
cols[is.na(cl)] <- "purple"

state1 <- c("1XCM_A", "3KKN_A", "4EFL_A", "4EFM_A", "4EFN_A")
cols[state1] <- "orange"

# some manul correction
cl["1XD2_A"] <- "TP"; cols["1XD2_A"] <- "red"
cl["1HE8_B"] <- "TP"; cols["1HE8_B"] <- "red"
cl["1WQ1_R"] <- "TP"; cols["1WQ1_R"] <- "red"
cl["4NMM_A"] <- "DP"; cols["4NMM_A"] <- "green"

## 1. 1XD2_A should be TP cause it has GDP+PO4
## 2. 1HE8_B should be TP (PDB annotation is wrong from pdb.annotate..)
## 3. 1WQ1_R change pink (transition state) to TP to make it simpler (cause only one pink..)
## 4. 4NMM_A the inhibitor ligand is similar to GDP (bind to the same site)
## 5. 4PZY_A also contains inhibitor 2XR
## 6. 5B30_A is reported as state1, but its structure similar to state2. I color it red now.

ras_annotate_trim <- cbind(ras_annotate_trim, "ligandType"=cl, "ligandColors"=cols)
write.csv(ras_annotate_trim, file="ras_annotate/ras_annotate_trim.csv")

# change data.frame into matrix
ras_annotate <- as.matrix(ras_annotate)
ras_annotate_all <- as.matrix(ras_annotate_all)
ras_annotate_trim <- as.matrix(ras_annotate_trim)

save(ras_annotate, ras_annotate_all, ras_annotate_trim,
     file="ras_annotate.RData")

#########
## pca ##
#########

# Find the subset of core residues
core <- core.find(pdbs_ras)
# Plot volume vs length
plot(core)
# Select core residues based on volume-length plot; default <= 1 A^3
core.inds <- print(core)

# Superpose the structures on the core residues
xyz_ras <- pdbfit(pdbs_ras, inds=core.inds, outpath="corefit/")


# Write core structure
write.pdb(xyz = xyz_ras[1, core.inds$xyz], resno=pdbs_ras$resno[1, core.inds$atom], file = "pca/ras_core.pdb")

# PCA on all structures
gaps.res <- gap.inspect(pdbs_ras$ali)
gaps.pos <- gap.inspect(pdbs_ras$xyz)

pca_ras <- pca.xyz(xyz_ras[, gaps.pos$f.inds])


# Write PC trajectory
a <- mktrj.pca(pca_ras, pc=1, file="pca/ras_pc1.pdb",
               resno = pdbs_ras$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_ras$ali[1, gaps.res$f.inds]) )
b <- mktrj.pca(pca_ras, pc=2, file="pca/ras_pc2.pdb",
               resno = pdbs_ras$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_ras$ali[1, gaps.res$f.inds]) )
c <- mktrj.pca(pca_ras, pc=3, file="pca/ras_pc3.pdb",
               resno = pdbs_ras$resno[1, gaps.res$f.inds],
               resid = aa123(pdbs_ras$ali[1, gaps.res$f.inds]) )

plot(pca_ras, col=cols[match(pdbs_ras_new$id, names(cols))])
dev.copy2pdf(file="figures/pca_ras.pdf")

save(core, core.inds, pca_ras, pdbs_ras, pdbs_ras_new, xyz_ras,
     file="pca.RData")




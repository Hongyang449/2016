## name: pdbs_annotation.r
## date: 09/12/2017

## Here I want to add more annotations to my pdbs

load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/ras_annotate.RData")

tbl_ras <- ras_annotate_trim[,c("ligandType","ligandColors")]

##############
## isoforms ##
##############

## detect hras/kras (no nras so far)
inds_h <- sort(c(grep("HR", ras_annotate_trim[,"compound"]), grep("H-R", ras_annotate_trim[,"compound"])))
inds_k <- sort(c(grep("KR", ras_annotate_trim[,"compound"]), grep("K-R", ras_annotate_trim[,"compound"])))
## 4 more pdbs have no ras isoform labels; but their sequence indicate they are hras
inds_h <- sort(c(inds_h, which(!(1:dim(ras_annotate_trim)[1] %in% c(inds_h,inds_k)))))
length(inds_h); length(inds_k)
# 115 hras; 28 kras

isoform <- rep(NA, dim(tbl_ras)[1])
isoform[inds_h] <- "H"
isoform[inds_k] <- "K"

tbl_ras <- cbind(tbl_ras, isoform)

##############
## mutation ##
##############

## summay of mutations ; ONLY for non-gap positions!
load("/Users/hyangl/project/ras/results/2016/0415_pdbs/ras/pdbs.RData")

hseq <- pdbs_ras_new$ali["5P21_A",!is.gap(pdbs_ras_new)]
kseq <- pdbs_ras_new$ali["3GFT_A",!is.gap(pdbs_ras_new)]
kseq[61] <- "Q"

hmut <- apply(pdbs_ras_new$ali[pdbs_ras_new$id[inds_h],!is.gap(pdbs_ras_new)], 1, function(x) {
  inds <- which(x != hseq)
  paste0(paste0(hseq[inds], inds, x[inds]), collapse=" ")
  })
kmut <- apply(pdbs_ras_new$ali[pdbs_ras_new$id[inds_k],!is.gap(pdbs_ras_new)], 1, function(x) {
  inds <- which(x != kseq)
  paste0(paste0(kseq[inds], inds, x[inds]), collapse=" ")
  })

mutation <- rep(NA, dim(tbl_ras)[1])
mutation[inds_h] <- hmut
mutation[inds_k] <- kmut

tbl_ras <- cbind(tbl_ras, mutation)

## distribution of mutations at each postion
mut_hras <- apply(pdbs_ras_new$ali[pdbs_ras_new$id[inds_h],!is.gap(pdbs_ras_new)], 2, table)
mut_kras <- apply(pdbs_ras_new$ali[pdbs_ras_new$id[inds_k],!is.gap(pdbs_ras_new)], 2, table)


#################
## extraLigand ##
#################

# inds for non-inhibitor structures
inds <- c(which(ras_annotate_trim[,"ligandId"] == "GTP,MG"),
         which(ras_annotate_trim[,"ligandId"] == "GDP,MG"),
         which(ras_annotate_trim[,"ligandId"] == "GNP,MG"),
         which(ras_annotate_trim[,"ligandId"] == "GCP,MG"),
         which(ras_annotate_trim[,"ligandId"] == "CA,GNP,MG"),
         which(ras_annotate_trim[,"ligandId"] == "PO4"),
         which(is.na(ras_annotate_trim[,"ligandId"])))

# manually check if ths structure contains an inhibitor
ras_annotate_trim[-inds,c("ligandId","ligandType")]

id_lig <- c("3L8Y_A", # different inhibitors
            "4DSO_A","4DST_A",
            "4EPT_A","4EPV_A","4EPW_A","4EPX_A","4EPY_A",
            "4LUC_A","4LUC_B","4LV6_A","4LV6_B",
            "4NMM_A",
            "4PZY_B","4PZZ_A",
            "5F2E_A",
            "3K8Y_A","3LBH_A","3LBI_A","3OIU_A","3OIW_A","4DLT_A","4DLU_A","4DLW_A", # CA ACT
            "3V4F_A","4DLR_A","4DLV_A","4DLX_A","4DLY_A","4DLZ_A") # DTT/DTU
# 2CE2, 2CL0, 2CL6, 2CL7, 2EVW are covalently linked to flurescent molucules; They are NOT treated as inhibitors here.
# different inhibitors; CA ACT; DTT/DTU
extraLigand <- rep("F", dim(tbl_ras)[1])

tbl_ras <- cbind(tbl_ras, extraLigand)
tbl_ras[id_lig, "extraLigand"] <- "T"

save(tbl_ras, mut_hras, mut_kras, inds_h, inds_k,
     pdbs_ras_new, ras_annotate_trim, 
     file="pdbs_annotation.RData")


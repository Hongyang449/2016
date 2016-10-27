## name: sse_gt.r
## date: 04/18/2016

attach(transducin)
sse <- bounds.sse(pdbs$sse[1,!is.gap(pdbs)])
# add length to sse for draw.box()
sse$helix$length <- sse$helix$end - sse$helix$start + 1
sse$sheet$length <- sse$sheet$end - sse$sheet$start + 1
# add sse to sse
sse$sse <- rep(" ",305)
for (i in 1:length(sse$helix$start)) {
  sse$sse[sse$helix$start[i]:sse$helix$end[i]] <- "H"
}
for (i in 1:length(sse$sheet$start)) {
  sse$sse[sse$sheet$start[i]:sse$sheet$end[i]] <- "E"
}

sse_gt <- sse

# sse for PL,SI,SII,a3
# use which(membership_gt==2)
helix <- list(start=6, end=13)
helix$length <- helix$end - helix$start + 1
PL <- list(helix=helix)
helix <- list(start=141, end=147)
helix$length <- helix$end - helix$start + 1
SI <- list(helix=helix)
helix <- list(start=164, end=182)
helix$length <- helix$end - helix$start + 1
SII <- list(helix=helix)
helix <- list(start=195, end=226)
helix$length <- helix$end - helix$start + 1
a3 <- list(helix=helix)

# add the concise one - excluding small unimportant sse
sse_gt_concise <- sse_gt
inds1 <- -c(3,5,7,10)
sse_gt_concise$helix$start <- sse_gt_concise$helix$start[inds1]
sse_gt_concise$helix$end <- sse_gt_concise$helix$end[inds1]
sse_gt_concise$helix$length <- sse_gt_concise$helix$length[inds1]
inds2 <- -c(1,6)
sse_gt_concise$sheet$start <- sse_gt_concise$sheet$start[inds2]
sse_gt_concise$sheet$end <- sse_gt_concise$sheet$end[inds2]
sse_gt_concise$sheet$length <- sse_gt_concise$sheet$length[inds2]


save(sse_gt, sse_gt_concise,
     PL, SI, SII, a3,
     file="sse_gt.RData")


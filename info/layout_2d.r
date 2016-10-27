## name: layout_2d.r
## data: 10/20/2016

load("/Users/hyangl/project/ras/results/2015/0730_ras_transducin/network/layout_2d.RData")
layout_ras <- layout_2d[1:9,]
layout_gt <- layout_2d
layout_eftu <- layout_2d; layout_eftu[11,] <- layout_eftu[11,] + c(5,3)

save(layout_ras, layout_gt, layout_eftu,
     file="layout_2d.RData")





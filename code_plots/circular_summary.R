library(circlize)
library(tidyverse)

# make a chrod diagram

# here you can check the similarity among mappes and callers

random_values <- c(500:100)
random_sample <- sample(random_values,12)
random_sample

col.pal = c(Freebayes="red",Lofreq="green",Varscan="blue",BreseQ="grey")

check <- matrix(random_sample,nrow=4,
       dimnames = list(c("Freebayes","Lofreq","Varscan", "BreseQ"),
                       c("breseq","BWA","Minimap")))

chordDiagram(check,grid.col = col.pal)



# make the plot more beautiful 
library(circlize)
circos.clear()
circos.initialize(factors=c("pb_501_001", "pb_501_002", "pb_501_003"), 
                  xlim=matrix(c(rep(0, 3), ref$V2), ncol=2))

# https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/

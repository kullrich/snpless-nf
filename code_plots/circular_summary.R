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

check

chordDiagram(check,grid.col = col.pal)


#https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/
  
library(tidyverse)
library(scales)
library(ggforce)

x1=seq(0,6000000,1000)
y1=rep(1.1,length(x1))
y2=rep(1.5,length(x1))

data2=as.data.frame(cbind(x1,y1,y2))
##
points= data.frame(mutations=c(0,5000000,10,1000000),check=c(1.7,1.7,1.3,1.3), type=c("indel","insertion"),
                   genome=c("leaky","leaky", "ANC", "ANC"), gene=c("wssE","wssA","pvdI","pvdII"))
points


ticks <- data.frame(
  x=seq(0,6000000,500000),
  y=rep(1,13),
  angle=360000-seq(0,360,30)
)



# organise and finsih plot to show it to Carsten
ggplot() +
  geom_rect(xmin=0, xmax=6000000, ymin= 1.1, ymax=1.5, fill="#eae3dc") +
  geom_line(data=data2,aes(x1,y1), color="#333333",size=1) +
  geom_line(data=data2,aes(x1,y2), color="#333333",size=1.5) +
  geom_point(data=points,aes(mutations,check, shape=type), size=5) +
  scale_shape_manual(values = c(15,18), name="Type of muatations") +
  coord_polar() +
  ylim(c(0,NA)) +
  annotate(geom="text", x=5500000, y=1.25, label="Ancestor", fontface=2, size=5, color="#666666") +
  annotate(geom="text", x=5500000, y=1.65, label="Leaky", fontface=2,size=5,color="#666666") +
  ggrepel::geom_label_repel(data=points, aes(mutations,check,label=gene),
                            box.padding = 0.5,
                            min.segment.length = 0,
                            nudge_x = .15,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  labs(title="Differences between colonies") +
  theme_void() +
  # theme(plot.title = element_text(hjust=0.5, face="bold", size = 20),
  #       panel.grid = element_blank(),
  #       axis.text.x= element_text(color="#333333", size=14),
  #       plot.background = element_rect(fill = "#eae3dc", color=  "#eae3dc"),
  #      panel.background = element_rect(fill ="#eae3dc", color=  "#eae3dc"), 
  #      axis.text.y=element_blank(),axis.ticks=element_blank(), 
  #      legend.background = element_rect(fill = "#eae3dc", color=  "#eae3dc"),
  #      legend.key  = element_rect(fill = "#eae3dc", color=  "#eae3dc"), 
  #      legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0,5000000,1000000), limits=c(0,6000000), 
                     labels = unit_format(unit = "MB", scale = 1e-6)) +
  annotate(geom="text", x=0, y=0, label="Pseudomonas fluorescens SBW25", 
           fontface="bold.italic",size=4,color="gray") +
  labs(x="",y="") +
  geom_text(data = ticks, aes(x, y, label = "|", angle = angle))

library('BioCircos')
myGenome = list("A" = 6000000)
check = list("B"=6000000)

BioCircos(genome = myGenome, yChr = FALSE, genomeFillColor = "Reds", chrPad = 0, 
          displayGenomeBorder = FALSE, genomeTicksDisplay = FALSE, genomeLabelDy = 0)

tracklist = BioCircosTextTrack('myTextTrack', 'Some text', size = "2em", opacity = 0.5, 
                               x = -0.67, y = -0.5)

BioCircos(tracklist, genomeFillColor = "PuOr",
          chrPad = 0, displayGenomeBorder = FALSE, 
          genomeTicksLen = 2, genomeTicksTextSize = 0, genomeTicksScale = 1e+8,
          genomeLabelTextSize = "9pt", genomeLabelDy = 0)

tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.5, maxRadius = 0.8,
                                     borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#FFBBBB")


BioCircos(tracklist, genomeFillColor = "PuOr",
          chrPad = 0.05, displayGenomeBorder = FALSE, 
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = "9pt", genomeLabelDy = 0)

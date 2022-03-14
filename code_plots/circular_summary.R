#__________________________
# ┬  ┬┌┐ ┬─┐┌─┐┬─┐┬┌─┐┌─┐
# │  │├┴┐├┬┘├─┤├┬┘│├┤ └─┐
# ┴─┘┴└─┘┴└─┴ ┴┴└─┴└─┘└─┘
#__________________________

## chord diagram
library(circlize)               # making circular plots
library(tidyverse)              # best package ever
library(ltc)                    # color palettes package


## plot genome differences 
library(scales)                 # you need to scale the x-y axis 
library(ggforce)                # to make shapes
library(ggalluvial)             # make an alluvial plot

#___________________________________________________________________________________
#
#  ██████╗██╗  ██╗ ██████╗ ██████╗ ██████╗     ██████╗ ██╗      ██████╗ ████████╗
# ██╔════╝██║  ██║██╔═══██╗██╔══██╗██╔══██╗    ██╔══██╗██║     ██╔═══██╗╚══██╔══╝
# ██║     ███████║██║   ██║██████╔╝██║  ██║    ██████╔╝██║     ██║   ██║   ██║   
# ██║     ██╔══██║██║   ██║██╔══██╗██║  ██║    ██╔═══╝ ██║     ██║   ██║   ██║   
# ╚██████╗██║  ██║╚██████╔╝██║  ██║██████╔╝    ██║     ███████╗╚██████╔╝   ██║   
#  ╚═════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═════╝     ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   
#___________________________________________________________________________________


# here you can check the similarity among mappers and callers
# create a test data frame

# ---------------------------
# create example data frame
# ---------------------------

test.df <- data.frame(Freebayes=sample(100:200,3), BCFtools=sample(100:200,3),
           Lofreq=sample(200:250,3), Varscan=sample(300:350,3), 
           Mappers=c("Breseq","BWA","Minimap")) %>% 
  pivot_longer(!Mappers, names_to = "SNP_caller")

# ---------------------------
# choose color -  palette
# ---------------------------

pltc(ltc("minou"))

# and choose the colors of your circus plot
col.pal = c(Freebayes="#00798c",Lofreq="#d1495b",Varscan="#edae49",BCFtools="#66a182",
            Breseq="#2e4057", BWA="#8d96a3", Minimap= "#c9cdd4")

# ---------------------------
# start making the plot
# ---------------------------

circos.clear()                # clear the area
                              # set the parameters of the chart
par(
  mar = c(1, 0, 3, 0),        # Margin around chart
  bg = c("white"),            # background color
  family="Avenir"
) 

                              # group the variables into callers and Mappers
SNP_callers = c("Freebayes","Lofreq","Varscan","BCFtools")
Mappers=c("Breseq","BWA","Minimap")
info=c(SNP_callers,Mappers)

# -----------------------------
# _|_ |_   _     _  _  ._ _  
#  |_ | | (/_   (_ (_) | (/_ 
# -----------------------------                 
                       
         
chordDiagram(test.df, order=info,
             directional = 1,            # ?
             #annotationTrack = "grid",   # ?
             diffHeight = F,            # ?
             preAllocateTracks = list(list(track.height=  uh(3,"mm")),   # outside track for names ?
                                      list(track.height=  uh(10,"mm"))), # middle track for regions ?
             self.link = 1,              # ?
             grid.col =col.pal,          # add colors in the bands
             transparency = 0.2)         # add tansparency in the bands 

# --------------------------------------                
# _|_ |_   _     _     _|_  _ o  _|  _    
#  |_ | | (/_   (_) |_| |_ _> | (_| (/_   
# --------------------------------------                                  
                                  
highlight.sector(SNP_callers, track.index = 1, col = "#ee6c4d",  # what is track index, I have no idea
                 text = "SNP_callers", cex = 0.7, text.col = "white", niceFacing = TRUE)

highlight.sector(Mappers, track.index = 1, col = "#3d5a80",
                 text = "Mappers", cex = 0.7, text.col = "white", niceFacing = TRUE)

# circos.track(track.index = 2, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=.7,col="grey75")
# }, bg.border = NA) # add state labels

# --------------------------------
# _|_ |_   _    _|_ o _|_ |  _    
#  |_ | | (/_    |_ |  |_ | (/_   
# --------------------------------                         
                          
title(main = list("Differences among callers and mappers",
                  cex=2.4,
                  col="#333333"))

#___________________________________________________________________________________________________________
# ╦═╗╔═╗╔═╗╔═╗╦═╗╔═╗╔╗╔╔═╗╔═╗╔═╗
# ╠╦╝║╣ ╠╣ ║╣ ╠╦╝║╣ ║║║║  ║╣ ╚═╗
# ╩╚═╚═╝╚  ╚═╝╩╚═╚═╝╝╚╝╚═╝╚═╝╚═╝

# 1. https://jokergoo.github.io/2020/05/21/make-circular-heatmaps/ [Make circular heatmaps with circlize]
# 2. https://r-charts.com/flow/chord-diagram/                      [Adjust colors in circlize ]
# 3. https://cran.r-project.org/web/packages/circlize/circlize.pdf [Manual of circlize, updates in June 2021]

#________________________________________________________________________________________________________
#
#  █████╗ ██╗     ██╗     ██╗   ██╗██╗   ██╗██╗ █████╗ ██╗         ██████╗ ██╗      ██████╗ ████████╗
# ██╔══██╗██║     ██║     ██║   ██║██║   ██║██║██╔══██╗██║         ██╔══██╗██║     ██╔═══██╗╚══██╔══╝
# ███████║██║     ██║     ██║   ██║██║   ██║██║███████║██║         ██████╔╝██║     ██║   ██║   ██║   
# ██╔══██║██║     ██║     ██║   ██║╚██╗ ██╔╝██║██╔══██║██║         ██╔═══╝ ██║     ██║   ██║   ██║   
# ██║  ██║███████╗███████╗╚██████╔╝ ╚████╔╝ ██║██║  ██║███████╗    ██║     ███████╗╚██████╔╝   ██║   
# ╚═╝  ╚═╝╚══════╝╚══════╝ ╚═════╝   ╚═══╝  ╚═╝╚═╝  ╚═╝╚══════╝    ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   
#________________________________________________________________________________________________________

pal=ltc("dora", 15, "continuous")

noPhage.colors <- function(unq.barcodes) {
  pal <- pnw_palette("Sailboat", n_distinct(test.df$value), type = "continuous")
  colorRampPalette(pal[1:n_distinct(test.df$value)])(length(levels(test.df$value)))
}

test.df

library("colorspace")

test.df %>% 
  #pivot_longer(cols = 1:2, names_to = "stratum", values_to = "program" ) %>% 
  rowid_to_column("alluvium") %>% 
  ggplot(aes(y=value, axis1=Mappers, axis2=SNP_caller)) +
  geom_alluvium(aes(fill=value)) +
  scale_fill_gradientn(colors=c(ltc("maya",12, "continuous")),
                       name="# of Variants",
                       guide=guide_legend(
                         direction = "horizontal",
                         title.position = "top"
                       )) +
  geom_stratum(width = 1/10, fill = "#2e4057", color = "black", ) +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)))
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), label.size=NA)   +
  scale_x_discrete(limits = c("Mappers", "SNP callers"), expand = c(.05, .05))  +
  ggtitle("Amount of SNPs that is called from combination of SNP callers and Mappers") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(family="Avenir", hjust=0.5),
        plot.title = element_text(hjust=0.5, family="Avenir", face="bold", size=18))

#___________________________________________________________________________________________________________
# ╦═╗╔═╗╔═╗╔═╗╦═╗╔═╗╔╗╔╔═╗╔═╗╔═╗
# ╠╦╝║╣ ╠╣ ║╣ ╠╦╝║╣ ║║║║  ║╣ ╚═╗
# ╩╚═╚═╝╚  ╚═╝╩╚═╚═╝╝╚╝╚═╝╚═╝╚═╝

# 1. https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html [Make alluvial plots with gg]
# 2. https://mstkolf.rbind.io/en/2020/08/02/tidytuesday-2020-week31/ [palmer penguins data]




#____________________________________________________________________________________________________
#
#  ██████╗██╗██████╗  ██████╗██╗   ██╗██╗      █████╗ ██████╗     ██████╗ ██╗      ██████╗ ████████╗
# ██╔════╝██║██╔══██╗██╔════╝██║   ██║██║     ██╔══██╗██╔══██╗    ██╔══██╗██║     ██╔═══██╗╚══██╔══╝
# ██║     ██║██████╔╝██║     ██║   ██║██║     ███████║██████╔╝    ██████╔╝██║     ██║   ██║   ██║   
# ██║     ██║██╔══██╗██║     ██║   ██║██║     ██╔══██║██╔══██╗    ██╔═══╝ ██║     ██║   ██║   ██║   
# ╚██████╗██║██║  ██║╚██████╗╚██████╔╝███████╗██║  ██║██║  ██║    ██║     ███████╗╚██████╔╝   ██║   
#  ╚═════╝╚═╝╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝    ╚═╝     ╚══════╝ ╚═════╝    ╚═╝   
#____________________________________________________________________________________________________

# ---------------------------
# create example data frame
# ---------------------------

genome_size = seq(0,6000000,1000)         # size of the genome
height1 = rep(1.1,length(genome_size))    # genome 1
height2 = rep(1.5,length(genome_size))    # genome 2

circles = as.data.frame(cbind(genome_size,height1,height2)) # compile it into a data.frame

variants1 = data.frame(var1.pos = sample(1:6000000,100), var1.height = 1.3)
variants2 = data.frame(var2.pos = sample(1:6000000,100), var2.height = 0.9)

ticks <- data.frame(    
  x1=seq(0,5500000,500000),
  y1=rep(0.7,12),
  y2=rep(0.6,12),
  x2=paste0((seq(0,5.5,0.5))),
  angle=360000-seq(0,330,30)
)

ggplot() +
  geom_rect(xmin=0, xmax=6000000, ymin= 1.1, ymax=1.5, fill="#eae3dc") +
  geom_line(data = circles , aes(genome_size,y1), color="#333333",size=3) +
  geom_line(data = circles , aes(genome_size,y2), color="#333333",size=3.5) +
  geom_jitter(data = variants1, aes(var1.pos, var1.height), size=3, color="#C02942", height=0.05, width = 0.05, alpha=0.7) +
  geom_jitter(data = variants2, aes(var2.pos, var2.height), size=3, color="#ECD078", height=0.05, width = 0.05, alpha=0.7) +
  coord_polar() +
  ylim(c(0,NA)) +
  geom_text(data = ticks, aes(x1, y1, label = "|", angle = angle)) +
  geom_text(data = ticks, aes(x1, y2, label = x2, angle = angle), size=5) +
  theme_void() +
  labs(title="Ploting variant positions across genome") +
  theme(plot.title = element_text(hjust=0.5))


point.mutations = data.frame(mutations=c(0,5000000,10,1000000), 
                   genome_coord=c(1.7,1.7,1.3,1.3),
                   mutation_type=c("indel","insertion"),
                   genome_name=c("leaky","leaky", "ANC", "ANC"), 
                   gene_name=c("wssE","wssA","pvdI","pvdII"))
point.mutations

#https://stackoverflow.com/questions/62556246/how-to-plot-the-variant-circular-bar-chart-in-r-with-ggplot

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


#_______________________________________________
# ┌─┐┌─┐┌┐┌┌─┐┌─┐┌┬┐┌─┐┌┐┌┌─┐┌┬┐┌─┐  ╦  ╦╔═╗╔═╗
# │  │ │││││  ├─┤ │ ├┤ │││├─┤ │ ├┤   ╚╗╔╝║  ╠╣ 
# └─┘└─┘┘└┘└─┘┴ ┴ ┴ └─┘┘└┘┴ ┴ ┴ └─┘   ╚╝ ╚═╝╚  
#_______________________________________________


library(here)
library(vcfR)
library(tidyverse)
library(strex)
library(data.table)
library(rvest)            # load rvest for reading html files
library(janitor)


setwd(here("code_plots"))
list.files()

vcf.raw=list.files(         # read in the directory all the html files
  pattern = "vcf", recursive = TRUE) 

vcf.merge = c()
for(i in 1:length(vcf.raw)) {
  vcf.dat1 <- read.vcfR(vcf.raw[i], verbose = FALSE) 
  id <- str_before_nth(vcf.raw[i], ".vcf", 1)        # make an id to keep only the first part of the names
  vcf.dat2 <- as.data.frame(vcf.dat1@fix) %>%    # you select from the VCF the fix option
    mutate(sample_id=id) %>% 
    mutate(POS=as.numeric(POS))
  
  vcf.merge <- rbind(vcf.merge,vcf.dat2)
} 

gd <- fread("BRESEQ.annotate.tsv")

html <- read_html("BRESEQ.annotate.html", skip=1) %>% 
  html_table(fill=TRUE, header = TRUE)  %>% 
  as.data.frame() %>% 
  janitor::row_to_names(row_number =1) %>% 
  filter(grepl("4,296,380",position))
html

  separate(.,position,into=c("position","note"),":", fill="right") %>% 
  mutate(position=as.numeric(str_replace_all(position,",",""))) 

vcf.merge %>% 
  filter(POS=="4296380")


html %>% 
  filter(position=="4296380")

tsv <- fread("BRESEQ.annotate.tsv") %>% 
  filter(position=="4296380") %>% 
  select(position,ref_seq,new_seq, title)
tsv






















# library('BioCircos')
# myGenome = list("A" = 6000000)
# check = list("B"=6000000)
# 
# BioCircos(genome = myGenome, yChr = FALSE, genomeFillColor = "Reds", chrPad = 0, 
#           displayGenomeBorder = FALSE, genomeTicksDisplay = FALSE, genomeLabelDy = 0)
# 
# tracklist = BioCircosTextTrack('myTextTrack', 'Some text', size = "2em", opacity = 0.5, 
#                                x = -0.67, y = -0.5)
# 
# BioCircos(tracklist, genomeFillColor = "PuOr",
#           chrPad = 0, displayGenomeBorder = FALSE, 
#           genomeTicksLen = 2, genomeTicksTextSize = 0, genomeTicksScale = 1e+8,
#           genomeLabelTextSize = "9pt", genomeLabelDy = 0)
# 
# tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0.5, maxRadius = 0.8,
#                                      borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#FFBBBB")
# 
# 
# BioCircos(tracklist, genomeFillColor = "PuOr",
#           chrPad = 0.05, displayGenomeBorder = FALSE, 
#           genomeTicksDisplay = FALSE,  genomeLabelTextSize = "9pt", genomeLabelDy = 0)

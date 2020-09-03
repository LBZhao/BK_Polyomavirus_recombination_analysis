# This script will create a circular diagram
# This script will take options
# example command in powershell: FOR %G IN (*.csv) DO Rscript "DIK_recombination_graph.R" %G
# -N  diagram of NCCR
# -T  diagram of total genome

library(ggraph)
library(igraph)
library(tidyverse)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
png(filename = paste(str_sub(args[1],end=-5),".png"), width = 1200, height = 1100, units = "px")

if (length(args)<2) {
  stop("Please provide Connection file and option -N for NCCR or -T for total genome", call.=FALSE)
Connection <- read.csv(file=args[1], header=FALSE, sep=",")
} else if (args[2] == "-N"){
d1=data.frame(from=c("origin", "origin", "origin", "origin", "origin"), to=c("O","P","Q", "R", "S"))
d2=data.frame(from="O", to=paste("subgroup", seq(1,142), sep="_"))
d3=data.frame(from="P", to=paste("subgroup", seq(143,210), sep="_"))
d4=data.frame(from="Q", to=paste("subgroup", seq(211,249), sep="_"))
d5=data.frame(from="R", to=paste("subgroup", seq(250,312), sep="_"))
d6=data.frame(from="S", to=paste("subgroup", seq(313,375), sep="_"))
hierarchy=rbind(d1, d2, d3, d4, d5, d6)
} else if (args[2] == "-T"){
  d1=data.frame(from=c("origin","origin","origin"), to=c("NCCR", "LATE", "EARLY"))
  d2=data.frame(from="NCCR", to=paste("subgroup", seq(1,375), sep="_"))
  d3=data.frame(from="LATE", to=paste("subgroup", seq(376,2676), sep="_"))
  d4=data.frame(from="EARLY", to=paste("subgroup", seq(2677,5141), sep="_"))
  hierarchy=rbind(d1, d2, d3, d4)
}

vertices = data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to))))


vertices$group = hierarchy$from[ match( vertices$name, hierarchy$to ) ]

mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

from = Connection$V1
to = Connection$V2

p=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  theme_void()

library(RColorBrewer)
if (args[2] == "-N"){
p +  geom_conn_bundle(data = get_con(from = from, to = to), width=0.9, alpha=0.02, colour="dodgerblue4", tension=0.9) +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   size=5) +
  scale_colour_manual(values= rep( brewer.pal(5,"Set1") , 30))
} else if (args[2] == "-T"){
p +  geom_conn_bundle(data = get_con(from = from, to = to), width=0.9, alpha=0.02, colour="dodgerblue4", tension=0.4) +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   size=5) +
  scale_colour_manual(values= rep( brewer.pal(3,"Set1") , 30))
}

dev.off()

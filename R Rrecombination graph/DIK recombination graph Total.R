# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
png(filename = paste(str_sub(args[1],end=-5),".png"), width = 1200, height = 1100, units = "px")
Connection <- read.csv(file=args[1], header=FALSE, sep=",")
#Connection <- data.frame(V1=c(5145),V2=c(5))
Connection$V3 <- with(Connection, V1+4)
Connection$V4 <- with(Connection, V2+4)

# create a data frame giving the hierarchical structure of your individuals
d1=data.frame(from=c("origin","origin","origin"), to=c("NCCR", "LATE", "EARLY"))
d2=data.frame(from="NCCR", to=paste("subgroup", seq(1,375), sep="_"))
d3=data.frame(from="LATE", to=paste("subgroup", seq(376,2676), sep="_"))
d4=data.frame(from="EARLY", to=paste("subgroup", seq(2677,5141), sep="_"))

hierarchy=rbind(d1, d2, d3, d4)

# create a vertices data.frame. One line per object of our hierarchy
vertices = data.frame(
  name = unique(c(as.character(hierarchy$from), as.character(hierarchy$to)))
) 

# Let's add a column with the group of each name. It will be useful later to color points
vertices$group = hierarchy$from[ match( vertices$name, hierarchy$to ) ]


# Create a graph object
mygraph <- graph_from_data_frame( hierarchy, vertices=vertices )

# The connection object must refer to the ids of the leaves:
from = Connection$V3
to = Connection$V4

p=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05)) +
  theme_void()

 #+  geom_conn_bundle(data = get_con(from = from, to = to), alpha=1, colour="RdPu", width=5, tension=1)



library(RColorBrewer)
p +  geom_conn_bundle(data = get_con(from = from, to = to), width=0.9, alpha=0.02, colour="dodgerblue4", tension=0.4) +
  geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group),   size=5) +
  scale_colour_manual(values= rep( brewer.pal(3,"Set1") , 30))
                                            
dev.off() 
                      
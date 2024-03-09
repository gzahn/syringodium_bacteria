# -----------------------------------------------------------------------------#
# Syringodium isoetifolium network analyses
# Exploring network connectivity
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     ggraph v 2.1.0
#                     ggmap v 3.0.1
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(ggraph); packageVersion("ggraph")
library(ggmap); packageVersion("ggmap")

# functions
source("./R/helper_functions.R")
source("./R/theme.R")
source("./R/googlemap_styling.R")

set.seed(666)

# Load google maps API key from .Renviron and set map style
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private
mapstyle <- rjson::fromJSON(file = "./R/mapstyle2.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)


# DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")
maxdist <- readRDS("./output/minimum_comm_dist_network.RDS")

# MAKE NETWORK ####
net <- make_network(ps %>% 
                    transform_sample_counts(function(x){x/sum(x)}),
                    keep.isolates = TRUE,max.dist = (1-maxdist),distance = "bray")

# PLOT NETWORK ####
plot_network(g = net,
             physeq = ps,
             color = "east_west",
             label = NULL,)

# simplified version
plot_net(ps %>% 
           transform_sample_counts(function(x){x/sum(x)}),
         color = "east_west",
         maxdist = maxdist,
         point_alpha = .5)
ggsave("./output/figs/simple_network_plot.png",dpi=300,height = 4,width = 4)

# merge by island before making network
ps@sam_data$merge_var <- paste(ps@sam_data$east_west,
                               ps@sam_data$location,
                               sep="_")
ps_island <- ps %>% 
  merge_samples("merge_var")
#repair metadata
ps_island@sam_data$east_west <- 
sample_names(ps_island) %>% 
  str_split("_") %>% 
  map_chr(1)
ps_island@sam_data$location <- 
sample_names(ps_island) %>% 
  str_split("_") %>% 
  map_chr(2)
sample_names(ps_island) <- ps_island@sam_data$location


net <- make_network(ps %>% 
                      transform_sample_counts(function(x){x/sum(x)}),
                    keep.isolates = TRUE,max.dist = 1-maxdist,distance = "bray")

# overlay on map
node.pos <- 
  data.frame(
  x = ps@sam_data$lon,
  y = ps@sam_data$lat
)
lay <- ggraph::create_layout(net,
                      layout = 'manual',
                      x=node.pos$x,
                      y=node.pos$y)
lay$name %>% str_split("_") %>% map_chr(1) %>% unique
# add weights and variables
igraph::V(net)
lay$weight <- igraph::degree(net,normalized = TRUE)
lay$east_west <- ps@sam_data$east_west

ge <- ggraph::get_edges()
lay_plot <- ge(lay)

lay_plot <- 
  lay_plot %>%
  mutate(relative_weight = round(node1.weight / max(node1.weight,na.rm = TRUE),2))


ggraph::ggraph(lay) +
  ggraph::geom_edge_link(aes(width=as.numeric(relative_weight/2)),data = lay_plot,width=1) +
  ggraph::geom_node_point(aes(color=east_west),size=3)

# PLOT NETWORK ON MAP ####
area <- 
  ggmap::get_googlemap(center = c(lon = ps@sam_data$lon %>% mean, 
                     lat = ps@sam_data$lat %>% mean),
        zoom = 5,
        scale = 2,
        style=mapstyle)

lay_plot$node2.name %>%  str_split("_") %>% map_chr(1) %>% unique

# map of sites
ggmap::ggmap(area) +
  geom_segment(aes(x=114.8,y=-11,xend=122,yend=6),linetype=2,linewidth=.25,color="red") +
  geom_point(data = microbiome::meta(ps_island), 
             aes(x=lon,y=lat,color=east_west),
             size=5,
             alpha=.75) +
  geom_text(data = microbiome::meta(ps_island), 
            aes(x=lon,y=lat,label=location),
            nudge_y = 1,
            size=2,
            fontface='bold') +
  scale_color_manual(values = pal.discrete,
                     breaks=c("West","East")) +
  labs(color="Side of\nWallace's Line",
       x="Longitude",
       y="Latitude") +
  theme(legend.position = "bottom",
        axis.title = element_text(face="bold",size=12))
ggsave("./output/figs/Location_Map.png",dpi=300,height = 6,width = 6)
ggsave("./output/figs/Figure_1.tiff",dpi=500,height = 6,width = 6)

wakatobi <- microbiome::meta(ps_island) %>% dplyr::filter(location=="Wakatobi")
banggai <- microbiome::meta(ps_island) %>% dplyr::filter(location=="Banggai")

lay_plot$weight <- lay_plot$relative_weight / nrow(lay_plot)

ggmap::ggmap(area) +
  geom_segment(aes(x=114.8,y=-11,xend=122,yend=6),linetype=2,linewidth=.25,color="red") +
  ggraph::geom_edge_link(data = lay_plot %>% mutate("Edge weight" = weight),
                         aes(colour=`Edge weight`),
                         show.legend=TRUE) +
  geom_point(data = lay_plot, aes(x=node1.x,
                                  y=node1.y,
                                  color = node1.east_west),
             size=5) +
  geom_point(data=wakatobi,aes(x=lon,y=lat),color=pal.discrete[2],size=5) + 
  geom_point(data=banggai,aes(x=lon,y=lat),color=pal.discrete[2],size=5) + 
  scale_color_manual(values = pal.discrete,
                     breaks=c("West","East")) +
  labs(color="Side of\nWallace's Line",
       x="Longitude",
       y="Latitude") +
  theme(legend.position = "right",
        axis.title = element_text(face="bold",size=12))

ggsave("./output/figs/mapped_network_plot_30maxdist.png",
       dpi=300,height = 6,width = 8)  

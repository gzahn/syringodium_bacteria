# Beta-Diversity
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(ggraph); packageVersion("ggraph")
library(ggmap); packageVersion("ggmap")
readRenviron("~/.Renviron")

source("./R/helper_functions.R")
source("./R/theme.R")
source("./R/googlemap_styling.R")

ps <- readRDS("./output/clean_phyloseq_object.RDS")

# Load google maps API key from .Renviron and set map style
ggmap::register_google(key = Sys.getenv("export APIKEY")) # Key kept private
mapstyle = 'feature:all|element:labels|visibility:off&style=feature:water|element:labels|visibility:on&style=feature:road|visibility:off'

# network analysis
net <- make_network(ps %>% 
                    transform_sample_counts(function(x){x/sum(x)}),
                    keep.isolates = TRUE,max.dist = .5,distance = "bray")

plot_network(g = net,
             physeq = ps,
             color = "east_west",
             label = NULL,)

# simplified version
plot_net(ps %>% 
           transform_sample_counts(function(x){x/sum(x)}),
         color = "east_west",
         maxdist = .7,
         point_alpha = .5)


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
                    keep.isolates = TRUE,max.dist = .3,distance = "bray")


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

# add weights and variables
lay$weight <- igraph::degree(net)
lay$east_west <- ps@sam_data$east_west
lay %>% head
ge <- ggraph::get_edges()
lay_plot <- ge(lay)

lay_plot <- 
  lay_plot %>%
  mutate(relative_weight = round(node1.weight / max(node1.weight,na.rm = TRUE),2))


ggraph::ggraph(lay) +
  ggraph::geom_edge_link(aes(width=as.numeric(relative_weight/2)),data = lay_plot,width=1) +
  ggraph::geom_node_point(aes(color=east_west),size=3)



mapstyle <- rjson::fromJSON(file = "./R/mapstyle2.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

area <- 
  ggmap::get_googlemap(center = c(lon = ps@sam_data$lon %>% mean, 
                     lat = ps@sam_data$lat %>% mean),
        zoom = 5,
        scale = 2,
        style=mapstyle)

lay_plot$weight <- lay_plot$relative_weight / nrow(lay_plot)
ggmap::ggmap(area) +
  geom_segment(aes(x=114.8,y=-11,xend=122,yend=6),linetype=2,linewidth=.25,color="red") +
  ggraph::geom_edge_link(data = lay_plot,
                         aes(alpha=weight),
                         show.legend=FALSE) +
  geom_point(data = lay_plot, aes(x=node1.x,
                                  y=node1.y,
                                  color = node1.east_west),
             size=5) +
  scale_color_manual(values = pal.discrete,
                     breaks=c("West","East")) +
  labs(color="Side of\nWallace's Line",
       x="Longitude",
       y="Latitude") +
  theme(legend.position = "bottom",
        axis.title = element_text(face="bold",size=12))

ggsave("./output/figs/mapped_network_plot_30maxdist.png",dpi=300,height = 6,width = 8)  

# -----------------------------------------------------------------------------#
# Syringodium isoetifolium distance decay
# Exploring distance decay of similarity (Mantel + MRM)
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     ecodist v 
#                     geodist v 
# -----------------------------------------------------------------------------#



# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(ecodist); packageVersion("ecodist")
library(geodist); packageVersion("geodist")

# DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS") %>% 
  transform_sample_counts(function(x){x/sum(x)})

# convert ASV table to data frame
otu <- ps@otu_table %>% as("matrix") %>% as.data.frame()

# pull sample data into data frame
meta <- ps@sam_data %>% as("data.frame")
class(meta)


# BRAY DISTANCE ####
# bray-curtis community distance matrix 
asv_dist <- vegdist(otu,method = "bray")

# SPATIAL DISTANCE ####
# lat-lon distance matrix
# nominal "distance" 
gps_dist <- vegdist(meta %>% select(lat,lon),
                    method = "euclidean")

# same thing, but in meters
haversine_dist <- geodist(meta %>% select(lon,lat),measure = "geodesic") %>% as.dist()

# MRM ####
# Multiple regression on matrices
mrm <- MRM(asv_dist ~ gps_dist, nperm = 1000)
mrm_haversine <- MRM(asv_dist ~ haversine_dist,nperm=1000)

# plot
data.frame(haversine = haversine_dist %>% as.matrix() %>% c(),
           asv = asv_dist %>% as.matrix() %>% c()) %>% 
  ggplot(aes(x=haversine,y=asv)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  labs(y="Community distance",x="Haversine distance (m)")
ggsave("./output/figs/comm_dist_vs_spatial_dist.png",dpi=300,height=4,width = 4)

# MANTEL ####
mant <- mantel(asv_dist,haversine_dist)
mant

# Mantel correlogram
correlog <- mantel.correlog(D.eco = asv_dist,
                D.geo = haversine_dist)
plot(correlog)



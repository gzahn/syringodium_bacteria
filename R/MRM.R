library(tidyverse)
library(ecodist)
library(phyloseq)
library(vegan)
library(geodist)

# load data
ps <- readRDS("./output/clean_phyloseq_object.RDS") %>% 
  transform_sample_counts(function(x){x/sum(x)})

# convert ASV table to data frame
otu <- ps@otu_table %>% as("matrix") %>% as.data.frame()

# pull sample data into data frame
meta <- ps@sam_data %>% as("data.frame")
class(meta)


# multipe regression on matrices ####
?MRM

# # examples
# data(graze)
# dist(graze$LOAR10)
# # Abundance of this grass is related to forest cover but not location
# MRM(dist(LOAR10) ~ dist(sitelocation) + dist(forestpct), data=graze, nperm=100)

# bray-curtis community distance matrix 
asv_dist <- vegdist(otu,method = "bray")

# lat-lon distance matrix
gps_dist <- vegdist(meta %>% select(lat,lon),
                    method = "euclidean")
# geodist (meters)
data.frame(lat = c(40.37463398553633,40.374667371814795),
           lon = c(-111.80117779863096,-111.80082938952233)) %>%
  geodist(measure = "geodesic") %>% 
  as.data.frame()

haversine_dist <- geodist(meta %>% select(lon,lat),measure = "geodesic") %>% as.dist()

# Multiple regression on matrices
mrm <- MRM(asv_dist ~ gps_dist, nperm = 1000)
mrm_haversine <- MRM(asv_dist ~ haversine_dist,nperm=1000)

# plot
data.frame(haversine = haversine_dist %>% as.matrix() %>% c(),
           asv = asv_dist %>% as.matrix() %>% c()) %>% 
  ggplot(aes(x=haversine,y=asv)) +
  geom_point() +
  geom_smooth(se=FALSE)

# Mantel test ####
mant <- mantel(asv_dist,haversine_dist)
mant

# Mantel correlogram
correlog <- mantel.correlog(D.eco = asv_dist,
                D.geo = haversine_dist)
plot(correlog)



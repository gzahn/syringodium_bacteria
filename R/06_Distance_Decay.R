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

set.seed(666)


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
mrm_haversine %>% as.data.frame() %>% 
  mutate(coef.asv_dist = coef.asv_dist %>% round(3),
         r.squared = r.squared %>% round(3)) %>% 
  saveRDS("./output/MRM_table.RDS")
# plot
data.frame(haversine = haversine_dist %>% as.matrix() %>% c(),
           asv = asv_dist %>% as.matrix() %>% c()) %>% 
  ggplot(aes(x=haversine,y=asv)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  labs(y="Community distance",x="Haversine distance (m)")
ggsave("./output/figs/comm_dist_vs_spatial_dist.png",dpi=300,height=4,width = 4)

# UNIFRAC DISTANCE
# GeoSphere sample distance vs UniFrac distance
lat <- ps@sam_data$lat
lon <- ps@sam_data$lon

# calculate positional distance between samples (relative, euclidean)
gps_dist <- 
  data.frame(lat,lon,row.names = sample_names(ps)) %>% 
  dist()

heatmap(as.matrix(gps_dist))

dist_m <- 
  geosphere::distm(data.frame(lon,lat)) %>% 
  as.dist()
UF <- UniFrac(ps %>% 
                transform_sample_counts(function(x){x/sum(x)}),
              weighted=TRUE)

data.frame(unifrac=c(UF),
           gps=c(dist_m)) %>% 
  ggplot(aes(x=dist_m,y=UF)) +
  geom_point(alpha=.05) +
  geom_smooth() +
  labs(x="Spatial distance (m)",
       y="Community distance (UniFrac)")

ggsave("./output/figs/comm_dist_vs_spatial_dist_unifrac.png",
       dpi=300,height = 6,width = 6)
# use GAM from MRM to set minimum community connectivity for network analysis
minimum_comm_dist_network <-   
  abs(
    data.frame(unifrac=c(UF),
               gps=c(dist_m)) %>%
      gam::gam(data = .,
               formula = unifrac ~ gps) %>% 
      coef() %>% 
      pluck("(Intercept)") - 1  
  )
saveRDS(minimum_comm_dist_network,
        "./output/minimum_comm_dist_network.RDS")


# MANTEL ####
mant <- vegan::mantel(asv_dist,haversine_dist)
call <- mant$call %>% as.character()

data.frame(
  call=paste0(call[1],"(xdis = ",call[2],", ydis = ",call[3],")"),
  test.stat=mant$statistic,
  P.value=mant$signif,
  N.perm=mant$permutations
  ) %>% 
  saveRDS("./output/mantel_table.RDS")

# Mantel correlogram
correlog <- mantel.correlog(D.eco = asv_dist,
                D.geo = haversine_dist)
plot(correlog)


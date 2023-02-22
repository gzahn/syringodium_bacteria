# -----------------------------------------------------------------------------#
# Syringodium isoetifolium beta diversity
# Exploring alpha diersity metrics
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
#                     patchwork v 1.1.2
#                     geodist v 0.0.8
#                     broom v 1.0.3
# -----------------------------------------------------------------------------#

# SETUP ####

# packages
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(geodist); packageVersion("geodist")
library(broom); packageVersion("broom")

# seed
set.seed(061223)

# data
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# ORDINATIONS ####

# Bray-Curtis NMDS
NMDS <- 
ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(method = "NMDS",
           distance = "bray")
plot_ordination(ps,NMDS,color="location") +
  stat_ellipse()

bray_plot <- plot_ordination(ps,NMDS,color="east_west") +
  stat_ellipse() +
  labs(title="Bray")

plot_ordination(ps,NMDS,color="lat")
plot_ordination(ps,NMDS,color="lon")



# Unifrac distance
# calculate weighted unifrac distance
UF <- UniFrac(ps %>% 
                transform_sample_counts(function(x){x/sum(x)}),
              weighted=TRUE)
# Do NMDS ordination on unifrac distance
UFORD <- ordinate(ps %>% 
           transform_sample_counts(function(x){x/sum(x)}),
         method = "NMDS",
         distance = UF)

# plot the unifrac ordination
unifrac_plot <- plot_ordination(ps,UFORD,color = "east_west") +
  stat_ellipse() +
  labs(title = "Unifrac")


bray_plot + unifrac_plot

# PERMANOVA ####

# specify model tables for Supplementary Info
vegan::adonis2(ps %>% 
                 transform_sample_counts(function(x){x/sum(x)}) %>% 
                 otu_table() ~ ps@sam_data$east_west + ps@sam_data$location) %>% broom::tidy() %>% 
  saveRDS("./output/bray_betadiv_df.RDS")
adonis2(UF ~ ps@sam_data$east_west + ps@sam_data$location) %>% broom::tidy() %>% 
  saveRDS("./output/unifrac_betadiv_df.RDS")

# beta-dispersion
bray <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  otu_table() %>% 
  vegdist()
vegan::betadisper(UF,ps@sam_data$east_west,type = 'centroid') %>% 
  saveRDS("./output/beta_dispersion_unifrac.RDS")



# send permanova results to text tables
sink("./output/permanova_tables.txt")
print("Bray distance")
vegan::adonis2(ps %>% 
                 transform_sample_counts(function(x){x/sum(x)}) %>% 
                 otu_table() ~ ps@sam_data$east_west)
print("Unifrac distance")
adonis2(UF ~ ps@sam_data$east_west + ps@sam_data$location)
sink(NULL)



# EXTRA STUFF? ####

lat <- ps@sam_data$lat
lon <- ps@sam_data$lon
sample_names(ps)

# calculate positional distance between samples (relative, euclidean)
gps_dist <- 
data.frame(lat,lon,row.names = sample_names(ps)) %>% 
  dist()

heatmap(as.matrix(gps_dist))

plot()

# plot distance matrices
data.frame(unifrac=c(UF),
           gps=c(gps_dist)) %>% 
  ggplot(aes(x=gps,y=unifrac)) +
  geom_point(alpha=.1) +
  geom_smooth(method="lm") +
  labs(title = "Unifrac distance vs positional distance")





# slope of linear model line
data.frame(unifrac=c(UF),
           gps=c(gps_dist)) %>% 
  glm(data = .,
      formula = unifrac ~ gps) %>% 
  coef()

dist_m <- 
geosphere::distm(data.frame(lon,lat)) %>% 
  as.dist()

data.frame(unifrac=c(UF),
           gps=c(dist_m)) %>% 
  ggplot(aes(x=dist_m,y=UF)) +
  geom_point() +
  geom_smooth(method="lm")
  

glm(data = .,
      formula = unifrac ~ gps) %>% 
  coef()


bk <- 
data.frame(lat,lon,location = sample_names(ps)) %>% 
  filter(grepl("Bali|Komodo",location))

geodist(bk %>% select(lon,lat)) %>% 
  dist()

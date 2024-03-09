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
library(lmerTest); packageVersion('lmerTest')

# seed
set.seed(061223)

# data
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# theme
theme_set(theme_minimal())
source("./R/theme.R")

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
  labs(title="Bray",color="Side of\nWallace's Line") +
  scale_color_manual(values = pal.discrete, breaks=c("West","East")) 

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
  labs(title = "Unifrac",color="Side of\nWallace's Line") +
  scale_color_manual(values = pal.discrete, breaks=c("West","East"))


bray_plot + unifrac_plot + 
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom') +
  theme(axis.title = element_text(face='bold',size=12),
        strip.text = element_text(face='bold',size=14),
        legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12))
ggsave("./output/figs/ordination_plots.png",dpi=300,height = 6, width = 12)
ggsave("./output/figs/Figure_3.tiff",dpi=500,height = 6, width = 12)

# PERMANOVA ####

# specify model tables for Supplementary Info

# Bray-Curtis permanova
vegan::adonis2(ps %>% 
                 transform_sample_counts(function(x){x/sum(x)}) %>% 
                 otu_table() ~ 
                 ps@sam_data$east_west + 
                 ps@sam_data$location,strata = ps@sam_data$east_west) %>% 
  broom::tidy() %>% 
  saveRDS("./output/bray_betadiv_df.RDS")

# UniFrac permanova
adonis2(UF ~ 
          ps@sam_data$east_west + 
          ps@sam_data$location,strata = ps@sam_data$east_west) %>% 
  broom::tidy() %>% 
  saveRDS("./output/unifrac_betadiv_df.RDS")

# beta-dispersion
bray <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  otu_table() %>% 
  vegdist()
vegan::betadisper(UF,ps@sam_data$east_west,type = 'centroid') %>% 
  saveRDS("./output/beta_dispersion_unifrac.RDS")



# # send permanova results to text tables
# sink("./output/permanova_tables.txt")
# print("Bray distance")
# vegan::adonis2(ps %>% 
#                  transform_sample_counts(function(x){x/sum(x)}) %>% 
#                  otu_table() ~ ps@sam_data$east_west)
# print("Unifrac distance")
# adonis2(UF ~ ps@sam_data$east_west + ps@sam_data$location)
# sink(NULL)



# EXTRA STUFF? ####



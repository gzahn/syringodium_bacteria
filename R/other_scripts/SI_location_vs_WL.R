
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


ps2 <- ps %>% 
  subset_samples(location %in% c("Bali","Tual","Pari","Derawan")) %>% 
  transform_sample_counts(function(x){x/sum(x)})

ord <- ps2 %>% 
  ordinate(method = "NMDS",
           distance = "bray")

plot_ordination(ps2,ord,color="location") +
  stat_ellipse() + theme(legend.position = 'none')



vegan::adonis2(ps2 %>% 
                 otu_table() ~ 
                 ps2@sam_data$east_west * 
                 ps2@sam_data$location) %>% 
  broom::tidy()

# -----------------------------------------------------------------------------#
# Syringodium isoetifolium differential abundance of taxa
# Exploring differential abundance of taxa
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.42.0
#                     Biostrings v 2.66.0
#                     corncob v 0.3.1
#                     vegan v 2.6.4
#                     patchwork v 1.1.2
#                     indicspecies v 1.7.12
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(indicspecies); packageVersion("indicspecies")


# Helper functions
source("./R/bbdml_helper.R")

# Random seed
set.seed(123)

# LOAD DATA ####
# Load cleaned phyloseq object
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# change order to west - east (since West feels more natural on left of plots)
ps@sam_data$east_west <- factor(ps@sam_data$east_west,levels = c("West","East"))

# merge taxa at genus level
ps_genus <- ps %>% 
  tax_glom(taxrank = "Genus")

# Clean up ASV names to show taxonomy
ASV_names <- otu_table(ps) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Phylum","Class","Order","Family","Genus"))

genus_names <- otu_table(ps_genus) %>% colnames()
genus_taxa <- otu_to_taxonomy(genus_names,ps_genus,level = c("Phylum","Class","Order","Family","Genus"))

# DIFFABUND TESTS ####

# use raw count data for corncob
da_analysis_eastwest <- differentialTest(formula = ~ east_west, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_genus,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
plot(da_analysis_eastwest) +
  theme(legend.position = 'none')

# pull model info out for reporting
mods <- da_analysis_eastwest$significant_models

capture_mods <- function(x){
 y <- x$coefficients %>% 
        as.data.frame() %>%
        janitor::clean_names()
 names(y) <- c("estimate","std_error","t_value","p_value")
 y <- y %>% 
   mutate(p_value = p_value %>% round(6))
 y <- y %>% 
   filter(row.names(y) %>% grepl(pattern="east_west"))
 return(y)
}
bbdml_mods <- map(mods,capture_mods)

# find the significant taxa
da_analysis_eastwest$significant_taxa
sig_taxa <- da_analysis_eastwest$significant_taxa %>% otu_to_taxonomy(data=ps_genus)

names(bbdml_mods) <- sig_taxa
joined_mods <- bbdml_mods %>% 
  purrr::reduce(full_join)
joined_mods$taxon <- names(bbdml_mods)
joined_mods <- joined_mods %>% 
  select(taxon,estimate,std_error,t_value,p_value)
joined_mods %>%
  saveRDS("./output/bbdml_significant_mod_tables.RDS")

# run bbdml() on all significant taxa
bbdml_obj <- multi_bbdml(da_analysis_eastwest,
                         ps_object = ps_genus,
                         mu_predictor = "east_west",
                         phi_predictor = "east_west",
                         taxlevels = 6)


# # filter significant taxa to only those with absolute effect sizes > 1.5
# is that unit log odds???

find_mu <- function(x,mu=1.5){
  y <- x$b.mu[2]
  z <- abs(y) > 1.5
  print(z)
}

new_bbdml_obj <- bbdml_obj[map(bbdml_obj, find_mu) %>% unlist]
names(new_bbdml_obj)


# INDICSPECIES ####
comm <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% as("matrix")
base::colnames(comm) <- genus_taxa

groups <- ps_genus@sam_data$east_west

# multipatt analysis
indval <- indicspecies::multipatt(comm,
                                  groups,
                                  control = how(nperm=999))

summary(indval)

# pull genus names of taxa that indicate EITHER East or West
x <- row.names(indval$A)
indicspecies_taxa <- x[which(indval$sign$p.value < 0.05)] %>% 
                      str_split("_") %>% 
                      map_chr(5)

# COMMON DIFFABUND TAXA ####
# find taxa that were detected by all methods
corncob_taxa <- sig_taxa %>% 
  str_split("_") %>% 
  map_chr(6)

# subset bbdml tests to just those found by multipatt as well
new_bbdml_obj <- new_bbdml_obj[which(names(new_bbdml_obj) %in% indicspecies_taxa)]
new_bbdml_obj %>% names() %>% 
  saveRDS("./output/final_significant_taxa.RDS")

# plot all significant taxa with mu > 1.5 and found by multipatt
plot_multi_bbdml(new_bbdml_obj,
                 color="east_west", 
                 pointsize = 3)

# save plots as RDS objects
plots <- ls(pattern = "bbdml_plot_")
for(i in plots){
  x <- base::get(i)
  saveRDS(x,file=file.path("./output/figs",paste0(i,".RDS")))
}

# reload objects (to avoid having to run plot_multi_bbdml() again)
# plotfiles <- list.files("./output/figs",pattern = "bbdml_plot_",full.names = TRUE)
# for(i in plotfiles){
#   bn <- base::basename(i) %>% str_remove(".RDS")
#   x <- readRDS(i)
#   assign(bn,x,envir = .GlobalEnv)
# }

# PLOT ####
# stick all plots together
plots <- ls(pattern = "^bbdml_plot_")

p1 <- bbdml_plot_1 +
  theme(legend.position = 'none',
        axis.title.y = element_blank())
p2 <- bbdml_plot_2 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p3 <- bbdml_plot_3 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') +
  labs(color="East or West of\nWallace's Line")
p4 <- bbdml_plot_4 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p5 <- bbdml_plot_5 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p6 <- bbdml_plot_6 +
  theme(legend.position = 'none',
        axis.title.y = element_blank()) +
  labs(color="East or West of\nWallace's Line")
p7 <- bbdml_plot_7 +
  theme(legend.position = 'none') +
  labs(y="Relative\nabundance")
p8 <- bbdml_plot_8 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p9 <- bbdml_plot_9 +
  theme(legend.position = 'none',
        axis.title.y = element_blank()) +
  labs(color="East or West of\nWallace's Line")
p10 <- bbdml_plot_10 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p11 <- bbdml_plot_11 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p12 <- bbdml_plot_12 +
  theme(legend.position = 'none',
        axis.title.y = element_blank()) +
  labs(color="East or West of\nWallace's Line")
p13 <- bbdml_plot_13 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p14 <- bbdml_plot_14 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p15 <- bbdml_plot_15 +
  theme(legend.position = 'none',
        axis.title.y = element_blank()) +
  labs(color="East or West of\nWallace's Line") 
p16 <- bbdml_plot_16 +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') 
p17 <- bbdml_plot_17 +
  theme(legend.position = 'bottom',
        axis.title.y = element_blank()) +
  labs(color="East or West of\nWallace's Line")
p18 <- bbdml_plot_18 +
  theme(legend.position = 'none',
        axis.title.y = element_blank()) +
  labs(color="East or West of\nWallace's Line")




(p1 + p2 + p3) / (p4 + p5 + p6) / (p7 + p8 + p9) / (p10 + p11 + p12) / (p13 + p14 + p15) / (p16 + p17 + p18)
ggsave("./output/figs/differential_abundance_sig_taxa.png", height = 8, width = 14,dpi=300)

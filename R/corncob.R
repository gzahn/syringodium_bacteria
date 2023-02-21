
# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")

source("./R/bbdml_helper.R")

# Load cleaned phyloseq object ####

ps <- readRDS("./output/clean_phyloseq_object.RDS")


# Differential abundance/dispersion tests ####

# Clean up ASV names to show taxonomy
ASV_names <- otu_table(ps) %>% colnames()
ASV_taxa <- otu_to_taxonomy(ASV_names,ps,level = c("Phylum","Class","Order","Family","Genus"))

ps_genus <- ps %>% 
  tax_glom(taxrank = "Genus")
ntaxa(ps_genus)
# use non-transformed data!
set.seed(123)
da_analysis_colcolor <- differentialTest(formula = ~ east_west, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_genus,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
plot(da_analysis_colcolor) + 
  theme()

# find the significant taxa
da_analysis_colcolor$significant_taxa
da_analysis_colcolor$significant_taxa %>% otu_to_taxonomy(data=ps_genus)

# This is a helper function I wrote. It's found in "scripts/bbdml_helper.R" 
bbdml_obj <- multi_bbdml(da_analysis_colcolor,
                         ps_object = ps_genus,
                         mu_predictor = "east_west",
                         phi_predictor = "east_west",
                         taxlevels = 6)

x <- bbdml_obj[[1]]
abs(x$b.mu[2]) > 1.5

find_mu <- function(x){
  y <- x$b.mu[2]
  z <- abs(y) > 1.5
  print(z)
}
map(bbdml_obj, find_mu) %>% unlist
new_bbdml_obj <- bbdml_obj[map(bbdml_obj, find_mu) %>% unlist]
new_bbdml_obj %>% names
# another helper function found in the same file
plot_multi_bbdml(bbdml_obj[1:3],
                 color="east_west", 
                 pointsize = 3)
plot_multi_bbdml
bbdml_plot_1

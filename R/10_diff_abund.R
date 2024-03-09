# SETUP ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(ranger)
library(vip)
library(corncob)
ra <- function(x){x/sum(x)}

# data ####
ps <- readRDS("./Output/clean_phyloseq_object.RDS")
sam <- sample_data(ps) %>% as("data.frame")
asv <- ps %>% transform_sample_counts(ra) %>% otu_table() %>% as('matrix')
asv <- as.data.frame(asv)
# colnames(asv) <- otu_to_taxonomy(taxa_names(ps),data = ps)
asv$sample_type <- sam$sample_type

dat <- asv %>% 
  mutate(sample_type = case_when(sample_type == "endophyte" ~ TRUE,
                                 TRUE ~ FALSE))
# ranger model
rf_mod <- ranger(data= dat,
                 formula = sample_type ~ .,importance = 'permutation')
imp <- vi(rf_mod)
imp$Variable %>% class
imp$Taxon <- otu_to_taxonomy(imp$Variable,data = ps)

imp %>% arrange(desc(Importance)) %>% head(10) %>% 
  ggplot(aes(x=Importance,y=Taxon)) +
  geom_col()

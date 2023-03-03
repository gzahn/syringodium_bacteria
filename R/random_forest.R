# Random Forest modeling
# can we predict side of wallaces line from taxa?

# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ranger); packageVersion("ranger")
library(ALEPlot); packageVersion("ALEPlot")
library(vip); packageVersion("vip")
library(corncob); packageVersion("corncob")

# DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS") # change to non-phylogeny stuff
ps_genus <- ps %>% 
  tax_glom("Genus") %>% 
  transform_sample_counts(function(x){x/sum(x)})

# get data ready for RF modeling
meta <- microbiome::meta(ps_genus) %>% 
  select(east_west) %>% 
  mutate(east = case_when(east_west == "East" ~ TRUE, # east is logical response
                          east_west == "West" ~ FALSE)) %>% 
  select(east)
asv <- otu_table(ps_genus) %>% 
  as("matrix") %>% 
  as.data.frame()

df <- 
meta %>% 
  bind_cols(asv)


# train RF model on full data set
ranger_model <- ranger::ranger(east~., data = df, classification = TRUE, probability = TRUE,importance = 'permutation')

# find top important taxa
top <- vip::vip(ranger_model) # find most important factors for success (survival)
top +
  theme(axis.text.y = element_blank())
vip_taxa <- corncob::otu_to_taxonomy(top$data$Variable,data = ps_genus) %>% unname()
vip_taxa

df[1,1]
pred <- predict(ranger_model,df) # how does it do predicting itself (no cross-validation)
preds <- pred$predictions %>% 
  as.data.frame() %>% 
  bind_cols(df$east)
names(preds) <- c("prob_T","prob_F","east")

preds %>% 
  ggplot(aes(x=prob_T,fill=east)) +
  geom_density()
vip_taxa %>% 
  saveRDS("./output/ranger_vip_taxa.RDS")

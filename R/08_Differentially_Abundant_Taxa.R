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
#                     microbiome v 1.20.0
#                     ranger v 0.14.1                    
#                     ALEPlot v 1.1
#                     vip v 0.3.2
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
library(microbiome); packageVersion("microbiome")
library(ranger); packageVersion("ranger")
library(ALEPlot); packageVersion("ALEPlot")
library(vip); packageVersion("vip")


# Helper functions
source("./R/bbdml_helper.R")
source("./R/theme.R")

# Random seed
set.seed(123)

# LOAD DATA ####
# Load cleaned phyloseq object
ps <- readRDS("./output/clean_phyloseq_object_silva.RDS")
ps@phy_tree <- NULL

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

# CORNCOB DIFFABUND ####

# At genus-level taxonomy (glommed)
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

# At ASV-level
# use raw count data for corncob
da_analysis_eastwest_ASV <- differentialTest(formula = ~ east_west, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
plot(da_analysis_eastwest_ASV) +
  theme(legend.position = 'none', axis.text.y = element_blank())

# pull model info out for reporting significant ASVs
mods <- da_analysis_eastwest_ASV$significant_models

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
da_analysis_eastwest_ASV$significant_taxa
sig_taxa <- da_analysis_eastwest_ASV$significant_taxa %>% otu_to_taxonomy(data=ps)

names(bbdml_mods) <- sig_taxa
joined_mods <- bbdml_mods %>% 
  purrr::map(as.data.frame) %>% 
  purrr::reduce(full_join)
joined_mods$taxon <- names(bbdml_mods)
joined_mods <- joined_mods %>% 
  select(taxon,estimate,std_error,t_value,p_value) %>% 
  mutate(ASV=da_analysis_eastwest_ASV$significant_taxa)
joined_mods %>%
  saveRDS("./output/bbdml_significant_mod_tables_ASV.RDS")


# export plot for Supplemental Info
plot(da_analysis_eastwest_ASV) +
  theme(legend.position = 'none', axis.text.y = element_blank())
ggsave("./output/figs/ASV-level_Diffabund_Plot.png",dpi=300,height=6,width=10)



# run bbdml() on all significant taxa
bbdml_obj <- multi_bbdml(da_analysis_eastwest,
                         ps_object = ps_genus,
                         mu_predictor = "east_west",
                         phi_predictor = "east_west",
                         taxlevels = 6)


# # filter significant taxa to only those with absolute effect sizes > 1.5 (log-odds)


find_mu <- function(x,mu=1.5){
  y <- x$b.mu[2]
  z <- abs(y) > 1.5
  print(z)
}

new_bbdml_obj <- bbdml_obj[map(bbdml_obj, find_mu) %>% unlist]
names(new_bbdml_obj)


# INDICSPECIES ####
comm <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% as("matrix")
base::colnames(comm) <- ASV_taxa

groups <- ps@sam_data$east_west

# deal with duplicated column names (multiple ASVs with same taxon names)
colnames(comm) <- paste(1:length(colnames(comm)),colnames(comm),sep=";")
# multipatt analysis
indval <- indicspecies::multipatt(comm,
                                  groups,
                                  control = how(nperm=999))

summary(indval)

# pull genus names of taxa that indicate EITHER East or West
x <- row.names(indval$A)


indicspecies_taxa <- x[which(indval$sign$p.value < 0.05)] %>%
                      str_replace("incertae_sedis","incertaesedis") %>% 
                      str_split("_") 
indicspecies_taxa <- lapply(indicspecies_taxa, `length<-`, max(lengths(indicspecies_taxa)))
indicspecies_taxa <- indicspecies_taxa %>% 
                       map_chr(5)
# RANDOM FOREST ####
# use relative abundance transformations for Ranger
ps_genus <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)})

# get data ready for RF modeling
meta <- microbiome::meta(ps) %>% 
  select(east_west) %>% 
  mutate(east = case_when(east_west == "East" ~ TRUE, # east is logical response
                          east_west == "West" ~ FALSE)) %>% 
  select(east)
asv <- otu_table(ps) %>% 
  as("matrix") %>% 
  as.data.frame()
df <- 
  meta %>% 
  bind_cols(asv)
# train RF model on full data set
ranger_model <- ranger::ranger(east~., 
                               data = df, 
                               classification = TRUE, 
                               probability = TRUE,
                               importance = 'permutation')

# how many top "vip" taxa should we select?
# find point at which variable importance from one place to the next drops (inflection point of exponential curve)
10^(lag(vip::vi(ranger_model)$Importance) - vip::vi(ranger_model)$Importance) %>% plot
(vip::vi(ranger_model)$Importance %>% cumsum() %>% round(3) -
    vip::vi(ranger_model)$Importance %>% cumsum() %>% round(3) %>% lag()) %>% head(30)

# top 20 capture all importance from the taxa (variables)
vip::vi(ranger_model) %>% 
  ggplot(aes(x=1:ntaxa(ps),y=Importance)) +
  geom_point() +
  geom_vline(xintercept = 20,linetype=2,color="red")

# find top important taxa with discriminatory power for model 
top <- vip::vip(ranger_model,num_features=20) # find most important factors for success (survival)
top +
  theme(axis.text.y = element_blank())
vip_taxa <- corncob::otu_to_taxonomy(top$data$Variable,data = ps_genus) %>% unname()
vip_taxa

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


# COMMON DIFFABUND TAXA ####
# find taxa that were detected by all methods
corncob_taxa <- sig_taxa

# subset bbdml tests to just those found by multipatt as well
new_bbdml_obj <- new_bbdml_obj[which(names(new_bbdml_obj) %in% indicspecies_taxa)]

# also subset bbdml tests to those found by Ranger as well
vip_taxa <- vip_taxa %>% 
  str_replace("incertae_sedis","incertaesedis") %>% 
  str_split("_") 
vip_taxa <- lapply(vip_taxa, `length<-`, max(lengths(vip_taxa)))
vip_taxa <- vip_taxa %>% 
  map_chr(6)
new_bbdml_obj <- new_bbdml_obj[which(names(new_bbdml_obj) %in% vip_taxa)]
# new_bbdml_obj <- new_bbdml_obj[which(names(new_bbdml_obj) != "Afipia")]
# names(new_bbdml_obj) <- names(new_bbdml_obj) %>% str_replace( "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium","Rhizobium")
new_bbdml_obj %>% names() %>% 
  saveRDS("./output/final_significant_taxa.RDS")
new_bbdml_obj %>% 
  saveRDS("./output/final_significant_bbdml_list.RDS")

# plot all significant taxa with mu > 1.5 and found by multipatt
plot_multi_bbdml(new_bbdml_obj,
                 color="east_west", 
                 pointsize = 3)
bbdml_plot_1;bbdml_plot_2;bbdml_plot_3;bbdml_plot_4;bbdml_plot_5;bbdml_plot_6;bbdml_plot_7
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
  labs(color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(legend.position = 'none',
        axis.title.y = element_blank())
p2 <- bbdml_plot_2 +
  labs(color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p3 <- bbdml_plot_3 +
  labs(color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p4 <- bbdml_plot_4 +
  labs(color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
p5 <- bbdml_plot_5 +
  labs(y="Relative\nabundance",
       color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(legend.position = 'none')
p6 <- bbdml_plot_6 +
  labs(color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') 
p7 <- bbdml_plot_7 +
  labs(color="East or West of\nWallace's Line") +
  scale_color_manual(values = pal.discrete) +
  theme(axis.title.y = element_blank(),
        legend.position = 'none')
# p8 <- bbdml_plot_8 +
#   labs(color="East or West of\nWallace's Line") +
#   theme(axis.title.y = element_blank(),
#         legend.position = 'none')
# p9 <- bbdml_plot_9 +
#   labs(color="East or West of\nWallace's Line") +
#   theme(legend.position = 'none',
#         axis.title.y = element_blank())
# p10 <- bbdml_plot_10 +
#   labs(color="East or West of\nWallace's Line") +
#   theme(axis.title.y = element_blank())



fullplot <- 
(p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + plot_spacer()) +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')
ggsave("./output/figs/differential_abundance_sig_taxa.png", height = 8, width = 10,dpi=300)
ggsave(plot = fullplot, 
       filename = "./output/figs/Figure_4.tiff", 
       height = 8, width = 10,dpi=500,device = "tiff",units = "in",bg = 'white')

# Abundance of diffabund taxa
diff_taxa_samples <- 
ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(Genus %in% c("Mucilaginibacter",
                           "Rhizobium",
                           "Sphingomonas",
                           "Elizabethkingia",
                           "Bradyrhizobium",
                           "Marixanthomonas",
                           "Azonexus"))
psmelt(diff_taxa_samples) %>% 
  dplyr::select(Genus,east_west,Abundance) %>% 
  mutate(east_west = factor(east_west,levels=c("West","East"))) %>% 
  ggplot(aes(x=Genus,y=Abundance,fill=east_west)) +
  scale_fill_manual(values = pal.discrete) +
  geom_boxplot() +
  facet_wrap(~east_west,scales='free') +
  coord_flip() +
  labs(y="Relative abundance",fill = "Side of\nWallace's Line") +
  theme(axis.text.y = element_text(face='bold.italic'),
        legend.position = 'none') 
ggsave("./output/figs/SI_Fig_relabund_of_diff_taxa.png",height = 4,width = 8)

psmelt(diff_taxa_samples) %>% 
  dplyr::select(Genus,east_west,Abundance) %>% 
  mutate(east_west = factor(east_west,levels=c("West","East"))) %>% 
  group_by(Genus,east_west) %>% 
  summarize(mean_relabund = mean(Abundance),
            max_relabund = max(Abundance))

microbiome::meta(ps) %>% 
  select(east_west) %>% 
  bind_cols(sample_sums(diff_taxa_samples)) %>% 
  ggplot(aes(x=east_west,y=diff_taxa_sample_sums)) +
  geom_boxplot()

microbiome::meta(ps) %>% 
  select(east_west) %>% 
  bind_cols(diff_taxa_sample_sums) %>%
  group_by(east_west) %>% 
  summarize(meansum=mean(diff_taxa_sample_sums))

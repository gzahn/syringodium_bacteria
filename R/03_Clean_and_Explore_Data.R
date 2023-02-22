# -----------------------------------------------------------------------------#
# Syringodium isoetifolium phyloseq cleanup
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.40.0
#                     ShortRead v 1.54.0
#                     Biostrings v 2.64.0
#                     adegenet v 2 .1 10
#                     readxl v 1.4.1
#                     janitor v 2.1.0
#                     microbiome v 1.20.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#        Remove non-bacteria, chloroplast and mitochondrial sequences           #
#                                                                               #
#################################################################################

# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(adegenet); packageVersion("adegenet")
library(readxl); packageVersion("readxl")
library(janitor); packageVersion("janitor")
library(microbiome); packageVersion("microbiome")

source("./R/googlemap_styling.R")
# seed
set.seed(789)

# DATA ####
ps <- readRDS("./output/ps_not-cleaned_w_tree.RDS") # change to non-phylogeny stuff

# EXPLORE ####
ps_nonbact <- subset_taxa(ps, Kingdom != "Bacteria")

# quick plot to look at kingdom-level taxonomy
# ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
#   plot_bar(fill="Kingdom")
# ggsave("./output/figs/16S_Kingdom-Level_Taxonomic_Proportions.png",dpi=300) # save figure for later use

# same plot, but non-bacteria, for sanity check
# ps_nonbact %>% 
#   transform_sample_counts(function(x){x/sum(x)}) %>%
#   plot_bar(fill="Kingdom")


# REMOVE NON-BACTERIA ####
ps <- subset_taxa(ps, Kingdom == "Bacteria")
tax <- tax_table(ps)

ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./output/16S_ASV_reference_sequences.RDS")


# double-check metadata

# Bangka has both east and west?
bad_bangka <- which(ps@sam_data$east_west == "West" & ps@sam_data$location == "Bangka")
ps@sam_data$east_west[bad_bangka] <- "East" # fix it


# ADD EXTERNAL DATA ####

## seagrass genotypes (geoff) ####
gt <- read_xlsx("./data/SI_GenoTypes_All.xlsx") %>% 
  clean_names()
# clean up sample names to match ps object
sites <- gt$name %>% str_split("(?<=[A-Za-z])(?=[0-9])")
gt$name <- paste0(
  sites %>% map_chr(1),
  "_",
  sites %>% map_chr(2)
) %>% 
  str_replace("_1$","_01") %>% # add leading zeros
  str_replace("_2$","_02") %>% 
  str_replace("_3$","_03") %>% 
  str_replace("_4$","_04") %>% 
  str_replace("_5$","_05") %>% 
  str_replace("_6$","_06") %>% 
  str_replace("_7$","_07") %>% 
  str_replace("_8$","_08") %>% 
  str_replace("_9$","_09") %>% 
  str_replace("Derwan_","Derawan_") # fix differential spelling
# subset to only those that match ps object
gt <- gt[which(gt$name %in% sample_names(ps)), ]
# make genind object
genind <- df2genind(X = select(gt,starts_with("s")),
          ncode = 3,
          ind.names = gt$name,
          ploidy = 1)
# pop clust analysis
clust <- find.clusters(genind,n.pca = 20,max.n.clust = 12,n.clust = 12)
# build metadata
pop_group_metadata <- 
microbiome::meta(ps) %>% 
  full_join(data.frame(sample = names(clust$grp),
                       pop_group = clust$grp))
# add to ps object
x <- sample_data(pop_group_metadata) # convert to sample_data object
sample_names(x) <- x$sample # make sample names match
ps@sam_data <- x # add to ps

# Save RDS object for cleaned up Phyloseq object
saveRDS(ps, file = "./output/clean_phyloseq_object.RDS")



# MESSING ABOUT ####

# does pop_group predict microbiome?
asv <- 
ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table()
vegan::adonis2(asv ~ ps@sam_data$east_west + ps@sam_data$location + ps@sam_data$pop_group)

# how does population structure look on a map?
ggmap::register_google(key = Sys.getenv("export APIKEY")) # Key kept private
# custom map style from JSON file
mapstyle <- rjson::fromJSON(file = "./R/mapstyle2.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

area2 <- 
  ggmap::get_googlemap(center = c(lon = ps@sam_data$lon %>% mean, 
                       lat = ps@sam_data$lat %>% mean),
          zoom = 5,
          scale = 2,
          style=mapstyle3)
ggmap::ggmap(area2) 


meta <- microbiome::meta(ps)


ggmap::ggmap(area2) +
  geom_jitter(data=meta,
              aes(x=lon,y=lat,color=pop_group),
              width = .5,
              height = .5,
              size=2) +
  scale_color_viridis_d() +
  geom_segment(aes(x=114.8,y=-11,xend=122,yend=6),linetype=2,linewidth=.25) +
  labs(color="Host\npopulation group",
       caption = paste0("Population groupings based on ",
                        gt %>% select(-name) %>% ncol(),
                        " microsatellite markers."))

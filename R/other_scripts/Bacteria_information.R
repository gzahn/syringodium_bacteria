#------------------------------------------------#
# Possible source of bacterial trait information #
# The dataset is full of NAs; looking at species #
# could provide accurate information             #
# Package Versions: R v 4.1.1                    #
#                   tidyverse v 1.3.2            #
#------------------------------------------------#
# Load Libraries ####
library(tidyverse) 

# Load in data ####
NCBI_species <- read_rds("./data/NCBI_Species.Rds")
NCBI_Traits <- read_rds("./data/NCBI_Traits.Rds")

genera <- c("Mucilaginibacter", "Rhizobium", "Sphingomonas", "Acinetobacter", "Elizabethkingia", 
"Bradyrhizobium",   "Aureitalea",       "Rhodospirillum",   "Arcobacter",       "Marixanthomonas", 
"Pelagicoccus",     "Azonexus",         "Thiogranum" ,      "Bizionia",         "Paracoccus",      
"Vibrionimonas",    "Wandonia",        "Bauldia",          "Litorimonas",      "Anderseniella")

Study_genera <- NCBI_Traits %>% 
  filter(genus %in% genera) %>% 
  select(c(species,genus, metabolism, pathways, sporulation, motility, range_tmp,
           range_salinity, cell_shape, d1_lo, d1_up, d2_lo, d2_up, optimum_ph))

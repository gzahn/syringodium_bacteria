# -----------------------------------------------------------------------------#
# Syringodium isoetifolium 16S DADA2 pipeline
# Reassigning taxonomy using Silva v 131.1
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     dada2 v 1.24.0
#                     decontam v 1.16.0
#                     phyloseq v 1.40.0
#                     Biostrings v 2.64.0
#                     patchwork v 1.1.1
#                     readxl v 1.4.1
#                     janitor::clean_names() v 2.1.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")
library(readxl); packageVersion("readxl")

# load data
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# convert asv table to matrix
seqtab.nochim <- otu_table(ps) %>% as('matrix')





# ASSIGN TAXONOMY ####

# Use RDP training set for 16S
taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/silva_nr99_v138.1_train_set.fa.gz", multithread = (parallel::detectCores() - 4))

# Save intermediate taxonomy file
saveRDS(taxa, file = "./output/silva_Taxonomy_from_dada2.RDS")

# add_species names
taxa <- addSpecies(taxa, "./taxonomy/silva_species_assignment_v138.1.fa.gz")

# Save completed taxonomy file
saveRDS(taxa, file = "./output/silva_Taxonomy_from_dada2_sp.RDS")

# load phy tree (not based on taxonomic assignment algorithm results)
fit <- readRDS("./output/16S_fit_treeNJ_muscle.RDS")


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(ps)
fit$tree$tip.label <- taxa_names(otu)
phy <- phy_tree(fit$tree)

sample_names(met)
sample_names(otu)
taxa_names(otu)

taxa_names(phy)

ps2 <- phyloseq(otu,met,tax,phy)

# compare assignments btwn databases
sum(ps@tax_table[,6] == ps2@tax_table[,6],na.rm=TRUE) / ntaxa(ps2)

df <- data.frame(rdp=ps@tax_table[,6],
                 silva=ps2@tax_table[,6])
colnames(df) <- c('rdp','silva')

df <- 
  df %>% 
  mutate(consensus = case_when(is.na(rdp) & is.na(silva) ~ NA,
                               is.na(rdp) & !is.na(silva) ~ silva,
                               !is.na(rdp) & is.na(silva) ~ rdp,
                               rdp == silva ~ silva,
                               !is.na(rdp) & !is.na(silva) & rdp != silva ~ silva))
consensus <- df$consensus
names(consensus) <- row.names(df)
ps2@tax_table[,6] <- consensus



# save output
saveRDS(ps2,"./output/clean_phyloseq_object_silva.RDS")



# -----------------------------------------------------------------------------#
# Syringodium isoetifolium phyloseq cleanup
# Cleaning up the phyloseq object
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.40.0
#                     ShortRead v 1.54.0
#                     Biostrings v 2.64.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
#        Remove non-bacteria, chloroplast and mitochondrial sequences           #
#                                                                               #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")

# load phyloseq object with phylogenetic tree ####

ps <- readRDS("./output/ps_not-cleaned_w_tree.RDS") # change to non-phylogeny stuff

# Find non-bacteria ####
ps_nonbact <- subset_taxa(ps, Kingdom != "Bacteria")

# quick plot to look at kingdom-level taxonomy
ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")
ggsave("./output/figs/16S_Kingdom-Level_Taxonomic_Proportions.png",dpi=300) # save figure for later use

# same plot, but non-bacteria, for sanity check
ps_nonbact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar(fill="Kingdom")


# REMOVE NON-BACTERIA, CHLOROPLASTS, MITOCHONDRIA, and empty samples/taxa ####
ps <- subset_taxa(ps, Kingdom == "Bacteria")
tax <- tax_table(ps)

ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./output/16S_ASV_reference_sequences.RDS")

# Save RDS object for cleaned up Phyloseq object
saveRDS(ps, file = "./output/clean_phyloseq_object.RDS")
ps

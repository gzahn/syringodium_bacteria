
# rarefaction and rarefying ####
library(vegan)

set.seed(2)
m <- matrix(sample(0:4,16,replace=TRUE),nrow=4,
            dimnames = list(sample=paste0("sam_",1:4),
                            species=paste0("sp_",1:4)))
rarecurve(m,label = FALSE) # rarefaction curve
mr <- rrarefy(m,sample=4) # rarefy to 5 "reads" per sample
# look at raw and rarefied matrices
m
mr

# indicator species ####
library(indicspecies)
# species that can 'indicate' certain environments
# look for species that indicate East vs West



# dissimilarity ####
# bray-curtis vs jaccard index
m
vegdist(m,method = 'jaccard',binary = TRUE)

# NULL models ####
m
# randomly reshuffle data N times
# get distribution of random results
# compare real data to random dist.
group <- c("A","A","B","B")
adonis2(mr ~ group) # PermANOVA

# R2 vs P-value ####
# p-value (how significant?)
# r2 (amount of variance explained by term)
# p-value = estimate of how surprised you should be
#   to see the result you did see, given that the NULL
#   is actually false

# sequencing machines and why reverse reads suck ####
# chemicals crap out eventually
# must read one direction at a time

# phred scores drop at end of read because:
#  seq can get out of sync. The longer your read, the more
#  chance to get out of sync

# OTU vs ASV ####
# OTU = operational taxonomic unit
#   invented to deal with seq errors and/or 
#   intraspecific variation
# ASV = exact sequence variant
#   depends on having NO sequencing error
#   can inflate diversity if intraspecific variation

# ITS1 vs ITS2 ####
# fungal barcodes
# depends on primers and taxa of interest
# use both!?

# full eukaryote meta-amplicon (18S)
# fungi (ITS1/2)
# bacteria (16S .... v3, v4, v5, v6)



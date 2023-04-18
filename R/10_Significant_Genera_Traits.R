# -----------------------------------------------------------------------------#
# Syringodium isoetifolium bacterial genus traits
# Testing phenotype enrichment in significant taxa 
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     phyloseq v 1.42.0
#                     patchwork v 1.1.2
#                     BacDive v 0.8.0
# -----------------------------------------------------------------------------#

# SETUP ####

# Packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(patchwork); packageVersion("patchwork")
library(BacDive); packageVersion("BacDive")

# Random seed
set.seed(8675309)

# initialize bacdive connection (add your userid and password)
# BacDive URL: https://bacdive.dsmz.de/
bacdive <- BacDive::open_bacdive(username = Sys.getenv("BACDIVE_USER"),
                                 password = Sys.getenv("BACDIVE_PW"))


# bacdive helper function
get_bacdive_morphology <- function(genus){
  x <- BacDive::request(object = bacdive,
                        query = genus,
                        search = "taxon")
  # find bacdive id
  y <- x$results
  
  if(length(y) == 0){
    df <- data.frame(genus=genus,
                     length=NA,
                     width=NA,
                     shape=NA)
    return(df)
  }
  
  # get info on that id
  z <- BacDive::fetch(bacdive,y)
  
  # get morphology
  m <- z$results %>% 
    map("Morphology")
  
  m <- m %>% map("cell morphology")
  # subset to non-null entries
  m <- m[c(which(m %>% map(notnull) %>% unlist) %>% names)]
  
  # test to see if morphology present
  # if morphology of type strain is empty, check taxonomic synonyms
  if(length(m) == 0){
    
    # find synonyms
    taxonomy <- z$results %>% 
      map("Name and taxonomic classification")
    syn <- taxonomy %>% 
      map("LPSN") %>% 
      map("synonyms") %>% 
      map("synonym") %>% 
      unlist() %>% 
      unique()
    
    # if not synonyms, return NAs
    if(is.null(syn)){
      df <- data.frame(genus=genus,
                       length=NA,
                       width=NA,
                       shape=NA)
      return(df)
    } else {
      # re-run bacdive request with synonyms
      x <- BacDive::request(object = bacdive,
                            query = syn,
                            search = "taxon")  
      
      if(length(x$results) == 0){
        df <- data.frame(genus=genus,
                         length=NA,
                         width=NA,
                         shape=NA)
        return(df)
      } else {
        z <- BacDive::fetch(bacdive,x$results)
      }
      
    }
  }
  
  # get morphology again (perhaps using synonyms)
  m <- z$results %>% 
    map("Morphology")
  m <- m %>% map("cell morphology")
  # subset to non-null entries
  m <- m[c(which(m %>% map(notnull) %>% unlist) %>% names)]
  
  # test to see if morphology present
  # if morphology still not present, export NAs
  
  if(length(m) == 0){
    df <- data.frame(genus=genus,
                     length=NA,
                     width=NA,
                     shape=NA)
    return(df)
  } else {
    # get morphology
    m <- z$results %>% 
      map("Morphology")
    m <- m %>% map("cell morphology")
    # subset to non-null entries
    m <- m[c(which(m %>% map(notnull) %>% unlist) %>% names)]
    
    # cell length
    length <- m %>% 
      map("cell length")
    if(which(length %>% map(notnull) %>% unlist()) %>% length() ==0){
      length <- NA
    } else {
      length <- length[which(length %>% map(notnull) %>% unlist())] %>% 
        unlist()
    }
    
    # cell width
    width <- m %>% 
      map("cell width")
    if(which(width %>% map(notnull) %>% unlist()) %>% length() ==0){
      width <- NA
    } else {
      width <- width[which(width %>% map(notnull) %>% unlist())] %>% 
        unlist()
    }
    
    # cell shape
    shape <- m %>% 
      map("cell shape")
    if(which(shape %>% map(notnull) %>% unlist()) %>% length() ==0){
      shape <- NA
    } else{
      shape <- shape[which(shape %>% map(notnull) %>% unlist())] %>% 
        unlist()
    }
    
    l <- length(length)
    w <- length(width)
    s <- length(shape)
    max_rows <- min(c(l,w,s))
    
    length <- length[1:max_rows]
    width <- width[1:max_rows]
    shape <- shape[1:max_rows]
    
    df <- data.frame(genus=genus,
                     length=length,
                     width=width,
                     shape=shape)
    return(df)
  }
  
}

# LOAD DATA ####
# Load cleaned phyloseq object
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# load list of significant taxa
sig_taxa <- readRDS("./output/final_significant_taxa.RDS")

# get list of unique genera
genus_list <- ps@tax_table[,6] %>% 
  table()
genus_list <- genus_list %>% 
  as.data.frame() %>% 
  pluck(".") %>% 
  levels() 

# BacDive API ####

# build large list of BacDive cell morphology
morph_list <- list()

for(i in genus_list){
  df <- get_bacdive_morphology(i)
  morph_list[[i]] <- df
}  
saveRDS(morph_list,"./output/genus_morphology_data.RDS")
morph_list <- readRDS("./output/genus_morphology_data.RDS")

names(morph_list)
# morph_list <- readRDS("./output/genus_morphology_data.RDS")

# reduce to single data frame
morphology <- purrr::reduce(morph_list,full_join)

# ANALYSIS OF BACDIVE DATA ####
# minimum dimensions
morphology$min_length <- morphology$length %>% 
  str_remove(" µm") %>% 
  str_split("-") %>% 
  map_chr(1) %>% 
  as.numeric()

morphology$min_width <- morphology$width %>% 
  str_remove(" µm") %>% 
  str_split("-") %>% 
  map_chr(1) %>% 
  as.numeric()


# add signifigance column
morphology <- morphology %>% 
  mutate(signifigant_taxa = case_when(genus %in% sig_taxa ~ TRUE,
                                      TRUE ~ FALSE),
         min_dimension = ifelse((min_length - min_width) > 0, 
                                min_width,
                                min_length),
         est_vol = min_length * min_width,
         shape = shape %>% str_remove("-shaped"),
         est_surface_area = case_when(shape == "coccus" ~ 4*pi*((min_dimension/2)^2),
                                  shape != "coccus" ~ (2*pi*((min_length/2)^2))))

morphology %>% 
  dplyr::filter(signifigant_taxa) %>% 
  group_by(genus) %>% 
  summarize(Avg_Min = mean(min_length,na.rm=TRUE))
  write_csv("./output/significant_taxa_traits.csv")

  
morphology %>% 
  dplyr::filter(signifigant_taxa) %>% 
  mutate(length=length %>% str_remove(" µm") %>% str_split("-") %>% map_chr(1) %>% as.numeric,
         width=width %>% str_remove(" µm") %>% str_split("-") %>% map_chr(1) %>% as.numeric) %>% 
  group_by(genus) %>% 
  summarize(mean_length=mean(length,na.rm=TRUE),
            mean_width=mean(width,na.rm=TRUE))
  
# distribution plots
morphology %>% 
  ggplot(aes(x=min_dimension, fill = signifigant_taxa)) +
  geom_density(alpha=.5) +
  labs(x="Minimum dimension (µm)",
       fill="Significant\ntaxa")
ggsave("./output/figs/minimum_cell_dimension_distribution.png",
       dpi=300,height = 6, width = 6)

# GLM model output for minimum cell dimension
# minimum dimension is important for Reynold's Number
mod <- glm(data=morphology,
           formula = min_dimension ~ signifigant_taxa)
saveRDS(mod,"./output/cell_dimension_glm.RDS")
sink("./output/cell_min_dimension_glm_summary.txt")
mod %>% summary()
sink(NULL)
report::report(mod)
morphology
# Chi-Square test for enrichment in cell shape for significant taxa

xsqtest <- table(morphology$shape, morphology$signifigant_taxa) %>% 
  chisq.test()
saveRDS(xsqtest,"./output/cell_shape_xsq.RDS")

sink("./output/cell_shape_chi-sq_test.txt")
print("Genus significance and shape")
table(morphology$shape, morphology$signifigant_taxa)
print("")
table(morphology$shape, morphology$signifigant_taxa) %>% 
  chisq.test()
print("No apparent enrichment in cell shape morphology in significant taxa")
sink(NULL)

morphology$est_surface_area
table(morphology$shape)


# surface area 
morphology %>% 
  ggplot(aes(x=est_surface_area,fill=signifigant_taxa)) +
  geom_density()

mod_surfacearea <- 
  glm(data=morphology,
      formula = est_surface_area ~ signifigant_taxa)
report::report(mod_surfacearea)
mod_surfacearea %>% summary
sink("./output/cell_surface_area_glm_summary.txt")
summary(mod_surfacearea)
sink(NULL)


# export morphology database
write_csv(morphology)




morphology %>% 
  filter(genus == "Vibrionimonas")

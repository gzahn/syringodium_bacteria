# SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggmap); packageVersion("ggmap")
source("./R/googlemap_styling.R")
options(scipen = 999)

# maybe this whole thing should be a shiny app
# user could pick genus of interest and see the corresponding relabund map 
# avoid calling google API though... pre-save images, and have app call the files.


# prepare for maps
ggmap::register_google(key = Sys.getenv("export APIKEY")) # Key kept private
# custom map style from JSON file
mapstyle <- rjson::fromJSON(file = "./R/mapstyle2.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

# seed
set.seed(789)

# DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS") # change to non-phylogeny stuff
ps_genus <- ps %>% 
  tax_glom("Genus") %>% 
  merge_samples("location") %>% 
  transform_sample_counts(function(x){x/sum(x)})
sig_mods <- readRDS("./output/bbdml_significant_mod_tables.RDS")
final_sig_taxa <- readRDS("./output/final_significant_taxa.RDS")


# repair metadata
ps_genus@sam_data$location <- sample_names(ps_genus)

# subset to just significant (corncob) taxa  
ps_genus_sig <- 
  ps_genus %>% 
  subset_taxa(tax_table(ps_genus)[,6] %in% final_sig_taxa)
# pull out metadata for convenience
meta <- microbiome::meta(ps_genus_sig)

# Base map
area <- 
  ggmap::get_googlemap(center = c(lon = ps@sam_data$lon %>% mean, 
                                  lat = ps@sam_data$lat %>% mean),
                       zoom = 5,
                       scale = 2,
                       style=mapstyle)


# for-loop to plot all maps of each significant genus (saved as files)
for(i in seq_along(final_sig_taxa)){
  
  g <- sig_mods$taxon[i] %>% str_split("_") %>% map_chr(6)
  abund <- ps_genus_sig@otu_table[,i] %>% as("numeric")
  
  p <- ggmap(area) +
    geom_segment(aes(x=114.8,y=-11,xend=122,yend=6),linetype=2,linewidth=.25) +
    geom_point(data=meta,
               aes(x=lon,
                   y=lat,
                   color=abund),
               size=3) +
    labs(title = g) +
    theme(legend.title = element_blank()) +
    scale_color_viridis_c()

  ggsave(plot=p,
         filename=paste0("./output/figs/abundance_maps/","abund_map_of_",g,".png"),
         dpi=300,
         height = 5,
         width = 5,
         device = "png")
}

# add relabund data to metadata df
for(i in seq_along(final_sig_taxa)){
  
  g <- final_sig_taxa[i]
  abund <- ps_genus_sig@otu_table[,i] %>% as("numeric")
  meta[g] <- abund
  
}
# Plot significant taxa against longitude (red bar for wallace's line)
# note: these are just the significant taxa detected by both corncob (mu > 1.5) and indicspecies
# merged by location!

meta <- 
meta %>% 
  select(c("sample","location","east_west","lat","lon","blank","seqtab_rows","pop_group"),
         all_of(final_sig_taxa))
meta %>% 
  pivot_longer(-c("sample","location","east_west","lat","lon","blank","seqtab_rows","pop_group"),
               names_to = "Genus",values_to = "RA") %>% 
  ggplot(aes(x=lon,y=RA)) +
  geom_vline(xintercept = 117.5,linetype=2,color='red',linewidth=.5) +
  geom_point(alpha=.5) +
  geom_smooth(color='black',se=FALSE,alpha=.5) +
  facet_wrap(~Genus,scales = 'free_y') +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold.italic'),
        axis.title = element_text(face='bold',size=14),
        axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
  labs(y="Relative abundance\n",
       x="\nLongitude",
       caption = "Note: Free Y-Axis Scales\nRed line marks approximate location of 'Wallace's Line'")
ggsave("./output/figs/sig_taxa_relative_abundance_vs_longitude.png",dpi=300,
       height = 6,width = 10)


# samples not merged by location:
# ps_genus <- ps %>% 
#   tax_glom("Genus") %>% 
#   transform_sample_counts(function(x){x/sum(x)})
# meta <- microbiome::meta(ps_genus)
# # add relabund data to metadata df
# for(i in seq_along(sig_mods$taxon)){
#   
#   g <- sig_mods$taxon[i] %>% str_split("_") %>% map_chr(6)
#   abund <- ps_genus_sig@otu_table[,i] %>% as("numeric")
#   meta[g] <- abund
#   
# }
# 
# 
# 
# meta <- 
#   meta %>% 
#   select(c("sample","location","east_west","lat","lon","blank","seqtab_rows","pop_group"),
#          all_of(final_sig_taxa))
# meta %>% 
#   pivot_longer(-c("sample","location","east_west","lat","lon","blank","seqtab_rows","pop_group"),
#                names_to = "Genus",values_to = "RA") %>% 
#   ggplot(aes(x=lon,y=RA)) +
#   geom_vline(xintercept = 117.5,linetype=2,color='red',linewidth=.5) +
#   geom_point(alpha=.5) +
#   geom_smooth(color='black',se=FALSE,alpha=.5) +
#   facet_wrap(~Genus,scales = 'free_y') +
#   theme_minimal() +
#   theme(strip.text = element_text(face = 'bold.italic'),
#         axis.title = element_text(face='bold',size=14),
#         axis.text.x = element_text(angle=90,hjust=1,vjust=.5)) +
#   labs(y="Relative abundance\n",
#        x="\nLongitude",
#        caption = "Note: Free Y-Axis Scales\nRed line marks approximate location of 'Wallace's Line'")
# ggsave("./output/figs/sig_taxa_relative_abundance_vs_longitude_not-merged.png",dpi=300,
#        height = 6,width = 10)


# merge samples based on longitude

# make new merge variable (combo of others)
ps@sam_data$newvar <- 
  paste(ps@sam_data$east_west,
      ps@sam_data$location,
      ps@sam_data$lon,
      sep="_")

# merge on that new variable
ps_lon <- ps %>% 
  merge_samples("newvar")

# repair metadata that was lost, using the newvar name
ps_lon@sam_data$east_west <- 
  sample_names(ps_lon) %>% 
  str_split("_") %>% 
  map_chr(1)
ps_lon@sam_data$location <- 
  sample_names(ps_lon) %>% 
  str_split("_") %>% 
  map_chr(2)
ps_lon@sam_data$lon <- 
  sample_names(ps_lon) %>% 
  str_split("_") %>% 
  map_chr(3)
ps_lon@sam_data

# "melt" the phyloseq object into a giant data frame, while transforming to relabund
melt_lon <- ps_lon %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt()

# re-order the samples by converting to a factor and setting levels from west to east
melt_lon <- 
  melt_lon %>% 
  mutate(Sample = factor(Sample,
                          levels = melt_lon %>% 
                            arrange(lon) %>% 
                            pluck("Sample") %>% unique()))

# plot stacked bar chart (keep labels off since they're too much)
melt_lon %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Phylum)) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_point(aes(y=0,color=east_west),size=5) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(face='bold',size=16),
        legend.title = element_text(face='bold',size=16)) +
  scale_fill_viridis_d() +
  scale_color_manual(values = pal.discrete,breaks=c("West","East")) +
  labs(x="Locations (West to East)",
       y="Relative abundance (phylum)",
       color="West or East of\nWallace's Line") +
  guides(fill = "none")

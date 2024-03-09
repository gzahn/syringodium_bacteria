# -----------------------------------------------------------------------------#
# Syringodium isoetifolium alpha diversity
# Exploring alpha diersity metrics
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     vegan v 2.6.4
#                     phyloseq v 1.42.0
#                     broom v 1.0.3
#                     lmerTest v 3.1.3
# -----------------------------------------------------------------------------#

# SETUP ####

# packages 
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(lmerTest); packageVersion("lmerTest")

# functions
source("./R/helper_functions.R")
source("./R/theme.R")

# options
options(scipen=999)

# data
ps <- readRDS("./output/clean_phyloseq_object_silva.RDS")

m <- microbiome::meta(ps)
we <- m %>% 
  arrange(lon) %>% 
  pluck("location") %>% unique()
# ALPHA-DIV ESTIMATES ####
alpha <- estimate_richness(ps) %>% 
  select(Observed,Shannon, Simpson)
# convert island to factor for plot ordering
ps@sam_data$location <- factor(ps@sam_data$location,levels = we)

plot_richness(ps,
              x="location",
              measures = c("Observed","Shannon"),
              color="east_west") +
  labs(color="Side of\nWallace's Line",
       x="Island") +
  theme(legend.title = element_text(hjust=.5,size=14),
        legend.position = "bottom",
        strip.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12)) +
  scale_color_manual(values=pal.discrete,breaks = c("West","East"))
ggsave("./output/figs/alpha_diversity_east_west_island.png",dpi=300,
       height = 4,width = 4)
ggsave("./output/figs/Figure_2.tiff",dpi=500,
       height = 4,width = 4)

# add metadata columns to alpha data frame
alpha$Sample <- row.names(alpha)
alpha$east_west <- ps@sam_data$east_west
alpha$island <- ps@sam_data$location

# MODELING ####

# richness
# mod_observed <- glm(data=alpha,
#                     formula = Observed ~ east_west * island)

m2 <- lmer(data=alpha,
     formula = Observed ~ east_west + (1|island))
saveRDS(m2,"./output/lmer_mod_alpha.RDS")


sink("./output/alpha_mod_observed_summary.txt")
summary(m2)
report::report(m2)
sink(NULL)
lme4::lmer

broom::tidy(mod_observed) %>% 
  filter(p.value<0.05)

# shannon
mod_shannon <- lmer(data=alpha,
                   formula = Shannon ~ east_west + (1|island))


sink("./output/alpha_mod_shannon_summary.txt")
summary(mod_shannon)
report::report(mod_shannon)
sink(NULL)


# Merge at location level ####

ps@sam_data$newvar <- 
paste(
  ps@sam_data$location,
  ps@sam_data$east_west,
  sep="_"
)


ps_island <- 
merge_samples(ps, "newvar", fun = sum)
ps_island@sam_data
# repair metadata using sample_names
x <- sample_names(ps_island) %>% str_split("_")

ps_island@sam_data$location <- 
x %>% 
  map_chr(1)

ps_island@sam_data$east_west <- 
  x %>% 
  map_chr(2)

ps_island@sam_data$east_west <- 
ps_island@sam_data$east_west %>% factor(levels=c("West","East"))

ps_island %>% 
  tax_glom(taxrank = "Phylum") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_heatmap(taxa.label = "Phylum",sample.label = "location") + 
  facet_wrap(~east_west,scales = "free_x")
ggsave("./output/figs/phylum_heatmap_east-west.png",dpi=300,width = 6, height = 6)

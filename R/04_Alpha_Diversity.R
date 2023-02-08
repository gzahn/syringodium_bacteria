##### Alpha Diversity #####

# Dependencies:
# tidyverse v 1.3.2
# vegan v 2.6.4
# phyloseq v 1.42.0
# broom v 1.0.3

# SETUP ####

# packages 
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")

# functions
source("./R/helper_functions.R")
source("./R/theme.R")

# data
ps <- readRDS("./output/clean_phyloseq_object.RDS")


# Alpha diversity measures ####
alpha <- estimate_richness(ps) %>% 
  select(Observed,Shannon, Simpson)
plot_richness(ps,
              x="location",
              measures = c("Observed","Shannon"),
              sortby = "Shannon",
              color="east_west") +
  labs(color="Side of\nWallace's Line",
       x="Island") +
  theme(legend.title = element_text(hjust=.5),
        legend.position = "bottom") +
  scale_color_manual(values=pal.discrete)
ggsave("./output/figs/alpha_diversity_east_west_island.png",dpi=300,
       height = 4,width = 4)

# add metadata columns to alpha data frame
alpha$Sample <- row.names(alpha)
alpha$east_west <- ps@sam_data$east_west
alpha$island <- ps@sam_data$location

# MODELING ####

# glm richness
mod_observed <- glm(data=alpha,
                    formula = Observed ~ east_west * island)

sink("./output/alpha_mod_observed_summary.txt")
summary(mod_observed)
sink(NULL)


broom::tidy(mod_observed) %>% 
  filter(p.value<0.05)

# glm shannon
mod_shannon <- glm(data=alpha,
                    formula = Shannon ~ east_west * island)
broom::tidy(mod_shannon) %>% 
  filter(p.value<0.05)

sink("./output/alpha_mod_shannon_summary.txt")
summary(mod_shannon)
sink(NULL)

report::report(mod_observed)





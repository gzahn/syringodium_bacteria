library(tidyverse)
library(phyloseq)
library(vegan)
library(broom)

#functions
source("./R/helper_functions.R")
source("./R/theme.R")

#options
options(scipen=999)

#data

ps <- readRDS("./output/clean_phyloseq_object.RDS")


# Alpha Diversity measures ###

ps
alpha <- estimate_richness(ps) %>% 
  select(Observed,Shannon,Simpson)
View(alpha)



#add sample id column
alpha$Sample <- row.names(alpha)

alpha$east_west <- ps@sam_data$east_west
alpha$island <- ps@sam_data$location

#repair meteadata using sample_names
sample_names(ps_island)

## creating a box plot

alpha %>% group_by(island) %>% ggplot(aes(x = island, y = Observed, fill = east_west))+
  geom_boxplot()+
  labs(fill = "Side Of \n Wallace's Line", x = "Island", y = "Observed Alpha Diversity")+
  theme(legend.title = (element_text(hjust = .5)), axis.text = element_text(angle = 20, hjust = 0.5))+
  scale_fill_manual(values=pal.discrete)

ggsave("./output/figs/Alpha_Diversity_IEW.jpg")

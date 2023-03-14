##### Alpha Diversity #####

# Dependencies:
# tidyverse v 1.3.2
# vegan v 2.5.7
# phyloseq v 1.36.0
# broom v 0.7.10

# SETUP ####

# packages 
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(phyloseq); packageVersion("phyloseq")
library(broom); packageVersion("broom")
library(janitor)

# functions
source("./R/helper_functions.R")
source("./R/theme.R")

# DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")
dates <- read.csv("./data/syringodium_dates.csv")
uv <- read.csv("./data/uv_index_islands.csv")

# CLEANING ####

# cleaning dates
dates <- clean_names(dates)
dates$sample <- dates$library_name
dates <- dates %>% select("collection_date", "sample")

# converting to dataframe
df <- ps@sam_data %>% 
  as("data.frame")

# merging dates to df
df <- left_join(df, dates, by = "sample")


# META DATA ####

# What is our meta data?
# we are using UV index
# we have uv index data for each site for every day of the month it was collected
# missing information for Bali
# https://power.larc.nasa.gov/data-access-viewer/



# UV INDEX ####
# averaging uv index
# merging uv average to df
uv_merge <- uv %>% select("location", "uv_month_avg")
uv_merge <- unique(uv_merge)
df <- left_join(df, uv_merge, by = "location")


# FINAL CSV ####
# need a sample column and UV column
# exporting to csv to turn in 

final <- df %>% select("sample", "uv_month_avg")
write.csv(final, file = "./output/morelli_meta.csv", row.names = FALSE)

# EXPLORING ####
uv %>% ggplot(aes(x = DY, y = ALLSKY_SFC_UV_INDEX)) +
  geom_point() +
  facet_wrap(~ location) 


uv %>% ggplot(aes(x = DY, y = ALLSKY_SFC_UV_INDEX, color = location)) +
  geom_line() +
  scale_color_manual(values=m_palette)

df %>% ggplot(aes(x = location, y = uv_month_avg, color = location)) +
  geom_point() +
  facet_wrap(~ east_west) +
  scale_color_manual(values=m_palette)


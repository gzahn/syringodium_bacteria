# load packages ####
library(tidyverse)
library(readxl)
source("./R/theme.R")


# read and clean metadata ####
meta <- read_xlsx("./data/SI_Indo_Metadata_CorrectGPS.xlsx") %>% 
  janitor::clean_names() %>% 
  rename("east_west" = "east_or_west_of_wallace_line") %>% 
  separate(gps,into = c("lat","lon"),sep=" ") %>% 
  mutate(lat=as.numeric(lat),
         lon=as.numeric(lon)) %>% 
  select(-barcode) %>% 
  mutate(blank = grepl("^Blank",sample))

# clean up remaining "N/A" values
meta[meta == "N/A"] <- NA

# plots ####
meta %>% 
  filter(blank == FALSE) %>% 
  ggplot(aes(x=lon,y=lat,color=location)) +
  geom_point(size=10) +
  scale_color_manual(values=pal.discrete)


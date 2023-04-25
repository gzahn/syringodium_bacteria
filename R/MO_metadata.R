#--------------------------------------------------------------------------------------------------------------------------------#
# pH and salinity metadata                                                                                                                    #
# From: https://www.ncei.noaa.gov/archive/archive-management-system/OAS/bin/prd/jquery/accession/download/162565                 #
# Citation:                                                                                                                      #
# Olsen, A., R. M. Key, S. van Heuven, S. K. Lauvset, A. Velo, X. Lin, C. Schirnick, A. Kozyr, T. Tanhua, M. Hoppema,            #
# S. Jutterström, R. Steinfeldt, E. Jeansson, M. Ishii, F. F. Pérez and T. Suzuki. The Global Ocean Data Analysis Project        #
# version 2 (GLODAPv2) - an internally consistent data product for the world ocean, Earth System Science Data, 8, 297-323, 2016. #
# doi: 10.5194/essd-8-297-2016                                                                                                   #
#--------------------------------------------------------------------------------------------------------------------------------#

# Setup ####
# libraries
library(tidyverse)
library(janitor)

# read in data
ph <- readRDS("data/MO_metadata_raw.RDS")
seagrass <- read_csv("data/syringodium_dates.csv") %>% 
  clean_names()

# manipulate the pH data to get the values we want to look at ####
# get column names
colnames(ph)
# select columns with data we might care about
ph <- ph %>% 
  select("year", "month", "latitude", "longitude", "phts25p0", "phtsinsitutp", "salinity")
# phts25p0 = ph at total scale, 25 degrees C, and surface pressure (0dbar)
# phtsinsitutp = pH at total scale, in situ temperature and pressure

# filter by year to match the years our seagrass samples were collected
ph <- ph %>% 
  filter(year == c("2010","2011"))

# This is a global study, so we need to filter by latitude and longitude to get 
# close to where our seagrass samples were collected
# separate seagrass lat and lon for the purpose of determining boundaries
lat_lon <- seagrass %>% 
  separate(lat_lon, c("lat", "north_south", "lon", "east_west"), sep = " ")
# North boundary
north <- lat_lon %>% 
  filter(north_south == "N") %>% 
  arrange(desc(lat))
View(north)
# North boundary is 3.92645
# South boundary
south <- lat_lon %>% 
  filter(north_south == "S") %>% 
  arrange(desc(lat))
View(south)
# South boundary is 8.68645
# I'll use -8.68645 when filtering the pH data because that is how that spreadsheet identifies latitudes 
# on the south side of the equator
# East and West boundaries
# Being on the eastern side of the prime meridian, there are no "W" values in the east_west column
# East boundary
lat_lon <- lat_lon %>% 
  arrange(desc(lon))
View(lat_lon)
# East boundary is 132.63785
# West boundary
lat_lon <- lat_lon %>% 
  arrange(lon)
View(lat_lon)
# West boundary is 106.61062
# filter pH data by lat and lon
# there are values that fall within our latitude boundaries, but not any within our longitude boundaries
# I expanded them a little, the range ends up being from 136 to 140. 
# Adjusting the smaller longitude doesn't provide more values until you get to negative numbers
ph <- ph %>% 
  filter(between(latitude, -8.68645, 3.92645)) %>% 
  filter(between(longitude, 100, 140))

# remove pH lat and lon columns
# they served their purpose of getting close to our sample locations, we don't necessarily need them anymore
ph <- ph %>% 
  select(!c("latitude","longitude"))

# adjust month formatting
ph$month <- paste("0", ph$month, sep="")

# combine year and month columns
ph <- ph %>% 
  unite(collection_date, year, month, sep = "-")

# join pH and seagrass data on date column
df <- full_join(seagrass,ph, by = "collection_date")

# write to csv
df %>% 
  write_csv("./data/MO_metadata.csv")




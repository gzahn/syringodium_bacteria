#---------------------------------------------------------------------------------#
# Precipitation metadata                                                          #
# From: https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/locations/FIPS:ID/detail #
# Finds the total daily precipitation of the collection date
# Package versions: R v 4.1.1                                                     #
#                   tidyverse v 1.3.2                                             #
#                   geodist v 0.0.8                                               #
#---------------------------------------------------------------------------------#               

# SETUP ####
library(tidyverse)
library(phyloseq)
library(geodist)


# Load collection date data and clean ####
dates <- read_csv("./data/syringodium_dates.csv") %>% 
  filter(!is.na(Collection_Date))

# Separating long and lat into two columns
dates <- 
  dates %>% 
  mutate(North = grepl(x = Lat_Lon, pattern = "N"),
         MultiplyLat = case_when(North == TRUE ~ 1,
                          North == FALSE ~ -1)) %>% 
  separate(col = Lat_Lon, into = c("Lat", "Lon"), sep = "N | S")

# Change to numeric and further cleanup
dates$Lat <- as.numeric(dates$Lat)
dates$Lat <- dates$Lat * dates$MultiplyLat
dates$Lon <- dates$Lon %>% str_remove(" E")
dates$Lon <- as.numeric(dates$Lon)

# Delete unneeded columns
dates$North <- NULL
dates$MultiplyLat <- NULL
dates$geo_loc_name <- NULL
dates$Host <- NULL


# Load precipitation data and clean ####
precip_data <- read_csv("./data/Precip_raw.csv") %>% 
  rename(Collection_Date = DATE) %>% 
  select(c(NAME, LATITUDE, LONGITUDE, Collection_Date, PRCP))

# Match the format of dates to the metadata dates
# Create seperate columns for month and year
precip_data$Collection_Date <- format(precip_data$Collection_Date, "%Y-%m")
precip_data <- precip_data %>% 
  separate(col = Collection_Date,into = c("Year", "Month"),sep = "-")

# Change N/As to 0s for precipitation
precip_data <- 
  precip_data %>% 
  mutate(PRCP = case_when(is.na(PRCP) ~ 0,
                          TRUE ~ PRCP))
# Remove any dates that are don't match any dates in the metadata
#x <- dates$Collection_Date
#precip_data <- 
#  precip_data %>% 
#  filter(Collection_Date %in% x) %>% unique()

# Get the sum for precipitation for grouped by station and collection date
final_precip_df <- 
  precip_data %>% 
  group_by(NAME, Year, Month) %>% 
  summarize(Precipitation = sum(PRCP),
            LATITUDE = LATITUDE,
            LONGITUDE = LONGITUDE) %>%  unique()


# Finding the closest station to the collection sites ####

# Create a new data frame for the long/lat for both data frames
sam_data_lon <- dates$Lon %>% unique()
sam_data_lat <- dates$Lat %>% unique()
sam_data_coord <- data.frame(sam_data_lon, sam_data_lat)
precip_data_lon <- final_precip_df$LONGITUDE %>% unique()
precip_data_lat <- final_precip_df$LATITUDE %>% unique()
precip_data_names <- final_precip_df$NAME %>% unique()

# Adding a value that was dropped by unique to its correct place
precip_data_lat <- append(x = precip_data_lat,values = -8.217, after = 4)
precip_data_coord <- data.frame(precip_data_lon, precip_data_lat,precip_data_names)

# Use Geodist to find closest stations to the sampling location
distances <- 
  geodist(x = precip_data_coord,
          y = sam_data_coord,
          measure = "geodesic") %>% 
  as.data.frame()

# Manually getting the closer points
sam_data_coord$NAME <- c("ALOR MALI KALABAHI, ID",
                            "BUBUNG, ID",
                            "MENADO SAM RATULAN, ID", 
                            "H AS HANANDJOEDDIN, ID", 
                            "HASANUDDIN, ID",
                            "KALIMARU, ID",
                            "GALELA GAMARMALAMU, ID",
                            "SEMARANG, ID", 
                            "MUHAMMAD SALAHUDDIN, ID", 
                            "TAREMPA, ID", 
                            "SOEKARNO HATTA INTERNATIONAL, ID",
                            "DENPASAR NGURAH RAI, ID",
                            "KAIMANA, ID",
                            "BAU BAU BETO AMBIRI, ID")


# Merging the two dataseets ####
final_precip_df <- 
  final_precip_df %>% 
  select(NAME, Year, Month, Precipitation)

merged_intermediate <- merge(final_precip_df, sam_data_coord)
cleanedish_meta <- merge(merged_intermediate, dates)
cleanedish_meta <- cleanedish_meta %>% 
  separate(col = Collection_Date, into = c("Year2", "Month2"), "-")

clean_meta <- 
  cleanedish_meta %>% 
  filter(Lat == sam_data_lat,
         Year == Year2,
         Month == Month2)

clean_meta <- 
  clean_meta %>% 
  select(Precipitation,`Library Name`) %>% 
  rename(sample = `Library Name`)


# 192 samples in the ps@sam_data
# Adding in sample[13:16] for each site
new_row <-  list(Precipitation = 0.04,  sample = "Derawan_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.04,  sample = "Derawan_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.04,  sample = "Derawan_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.04,  sample = "Derawan_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Pari_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Pari_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Pari_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Pari_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Halmahera_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Halmahera_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Halmahera_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Halmahera_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Belitung_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Belitung_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Belitung_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Belitung_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Banggai_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Banggai_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Banggai_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Banggai_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Alor_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Alor_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Alor_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Alor_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Tual_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Tual_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Tual_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Tual_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Komodo_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Komodo_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Komodo_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Komodo_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Wakatobi_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Wakatobi_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Wakatobi_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Wakatobi_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 0.00,  sample = "Karimunjawa_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Karimunjawa_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Karimunjawa_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 0.00,  sample = "Karimunjawa_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = 1.54,  sample = "Bangka_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 1.54,  sample = "Bangka_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 1.54,  sample = "Bangka_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = 1.54,  sample = "Bangka_16")
clean_meta <-  rbind(clean_meta,new_row)

new_row <-  list(Precipitation = NA,  sample = "Bali_01")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_02")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_03")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_04")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_05")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_06")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_07")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_08")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_09")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_10")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_11")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_12")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_13")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_14")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_15")
clean_meta <-  rbind(clean_meta,new_row)
new_row <-  list(Precipitation = NA,  sample = "Bali_16")
clean_meta <-  rbind(clean_meta,new_row)

# Rows to filter out "bira", "natuna", "Sanur"
clean_meta <- clean_meta[!grepl("Bira", clean_meta$sample),]
clean_meta <- clean_meta[!grepl("Natuna", clean_meta$sample),]
clean_meta <- clean_meta[!grepl("Sanur", clean_meta$sample),]

# If you want to save the clean dataset
##write_csv(x = clean_meta, file = "./output/precip_metadata.csv")


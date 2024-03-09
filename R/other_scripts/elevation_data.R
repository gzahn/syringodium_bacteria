library(tidyverse)
library(googleway)

#Google API KEy
#Insert your google maps API here
set_key("")

#Read in Data

df <- read.csv("./data/syringodium_dates.csv", na.strings = c("","NA"))

#Remove NA rows

df <- df %>% filter(!is.na(Collection_Date))

#SPlitting latitude and longitude into two columns

x <- as.data.frame(str_split(df$Lat_Lon, " "),row.names = c("Lat","NS","Long","EW")) %>% t()

df1 <- as.data.frame(x, row.names = c(1:168))

df1 <- df1 %>% mutate(Lat=as.numeric(Lat), Long=as.numeric(Long))

#Converting NS lat/long format in +/- format

df1$NS <- case_when(df1$NS == "S" ~ -1, df1$NS == "N" ~ 1)

df1 %>% mutate(NS-as.numeric(NS))

df1$Lat <- df1$Lat * df1$NS

df1 <- df1 %>% select(c("Lat","Long"))

df$lat <- df1$Lat
df$lon <- df1$Long

#Calls Google earth API for elevation data on df


edata <- googleway::google_elevation(df_locations = df, location_type = "individual")
  


# Pulling out the elevation

df$elevation <- edata$results$elevation

# Writing to csv a dataframe with only sample name and 

elevation <- df %>% select("sample_name" = Library.Name,elevation)

elevation %>% write_csv("./Data/elevation_data.csv")



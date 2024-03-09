###############################################################
## Finding average air temperature per month at each station ##
## Author: Porter Bischoff                                   ##
###############################################################


# read in data
data2011 <- read.csv("./data/mean_air_temp_data/May 1 - Aug 31 2011.csv")
data2010 <- read.csv("./data/mean_air_temp_data/May 1 - Aug 31 2010.csv")
data_coll <- read.csv("./data/mean_air_temp_data/syringodium_dates.csv")

# load library ####
library(janitor)
library(tidyverse)
library(tidyr)
library(sf)
library(mapview)
library(dplyr)

# cleaning data2010 ####
# run through janitor
data2010 <- data2010 %>% clean_names()

# splitting date to different columns
data2010$date <- as.character(data2010$date)
data2010 <- separate(data2010, col=date, into=c('year', 'month', 'day'), sep='-')
data2010 <- separate(data2010, col=day, into = c('trash', 'time'), sep = 'T')
data2010 <- separate(data2010, col = time, into = c('hour','minute','second'), sep = ':')
data2010 <- subset(data2010, select = -c(trash, second, hour, minute))
data2010 <- data2010 %>% drop_na(air_temp)

#get station means
by_stat_10 <- data2010 %>% group_by(station)
air_temp_stat_10 <- by_stat_10 %>% 
  group_by(station) %>% 
  mutate(mean_station_temp = mean(air_temp))


# get monthly station means
by_month_10 <- by_stat_10 %>% group_by(station,month)
air_temp_month_10 <- by_month_10 %>% 
  group_by(station) %>% 
  group_by(month) %>% 
  mutate(mean_station_temp_month = mean(air_temp))

# cleaning data2011 ####
# run through janitor
data2011 <- data2011 %>% clean_names()

# splitting date to different columns
data2011$date <- as.character(data2011$date)
data2011 <- separate(data2011, col=date, into=c('year', 'month', 'day'), sep='-')
data2011 <- separate(data2011, col=day, into = c('trash', 'time'), sep = 'T')
data2011 <- separate(data2011, col = time, into = c('hour','minute','second'), sep = ':')
data2011 <- subset(data2011, select = -c(trash, second, hour, minute))
data2011 <- data2011 %>% drop_na(air_temp)

#get station means
by_stat_11 <- data2011 %>% group_by(station)
air_temp_stat_11 <- by_stat_11 %>% 
  group_by(station) %>% 
  mutate(mean_station_temp = mean(air_temp))


# get monthly station means
by_month_11 <- by_stat_11 %>% group_by(station,month)
air_temp_month_11 <- by_month_11 %>% 
  group_by(station) %>% 
  group_by(month) %>% 
  mutate(mean_station_temp_month = mean(air_temp))



# cleaning data_coll####
data_coll <- data_coll %>% clean_names()

# get rid of na values
data_coll[data_coll == "" | data_coll == " "] <- NA
df_coll <- na.omit(data_coll)

# splitting data to different columns
df_coll <- separate(df_coll, col = collection_date, into = c('year', 'month'), sep='-')
df_coll <- separate(df_coll, col = host, into = c('genus','species'), sep = ' ')
df_coll <- separate(df_coll, col = geo_loc_name, into = c('location','island'), sep = ': ')
df_coll <- subset(df_coll, select = -c(location))
df_coll <- separate(df_coll, col = lat_lon, into = c('lat','n_s','lon','e_w'), sep = ' ')
df_coll$lat <- as.numeric(df_coll$lat)
df_coll$lon <- as.numeric(df_coll$lon)


# make latitude negative at 'S'
df_coll <- df_coll %>% 
  mutate(latitude = case_when(
    n_s=='S' ~ lat*(-1),
    n_s=='N' ~ lat
  ))

df_coll <- df_coll %>% 
  mutate(longitude = case_when(
    e_w=='E' ~ lon,
    e_w=='W' ~ lon*(-1)
  ))
df_coll <- subset(df_coll, select = -c(lat, n_s, lon, e_w))

# separate the data into sections
df_coll <- df_coll %>% 
  mutate(section = case_when(
    longitude>=90 & longitude<105 & latitude>=-15 & latitude<(-8.75) ~ 1,
    longitude>=90 & longitude<105 & latitude>=-8.75 & latitude<(-2.5) ~ 2,
    longitude>=90 & longitude<105 & latitude>=-2.5 & latitude<3.75 ~ 3,
    longitude>=90 & longitude<105 & latitude>=3.75 & latitude<=10 ~ 4,
    longitude>=105 & longitude<120 & latitude>=-15 & latitude<(-8.75) ~ 5,
    longitude>=105 & longitude<120 & latitude>=-8.75 & latitude<(-2.5) ~ 6,
    longitude>=105 & longitude<120 & latitude>=-2.5 & latitude<3.75 ~ 7,
    longitude>=105 & longitude<120 & latitude>=3.75 & latitude<=10 ~ 8,
    longitude>=120 & longitude<135 & latitude>=-15 & latitude<(-8.75) ~ 9,
    longitude>=120 & longitude<135 & latitude>=-8.75 & latitude<(-2.5) ~ 10,
    longitude>=120 & longitude<135 & latitude>=-2.5 & latitude<3.75 ~ 11,
    longitude>=120 & longitude<135 & latitude>=3.75 & latitude<=10 ~ 12,
    longitude>=135 & longitude<=150 & latitude>=-15 & latitude<(-8.75) ~ 13,
    longitude>=135 & longitude<=150 & latitude>=-8.75 & latitude<(-2.5) ~ 14,
    longitude>=135 & longitude<=150 & latitude>=-2.5 & latitude<3.75 ~ 15,
    longitude>=135 & longitude<=150 & latitude>=3.75 & latitude<=10 ~ 16,
  ))

# separate into 'sections' ####
air_temp_month_10 <- air_temp_month_10 %>% 
  mutate(section = case_when(
    longitude>=90 & longitude<105 & latitude>=-15 & latitude<(-8.75) ~ 1,
    longitude>=90 & longitude<105 & latitude>=-8.75 & latitude<(-2.5) ~ 2,
    longitude>=90 & longitude<105 & latitude>=-2.5 & latitude<3.75 ~ 3,
    longitude>=90 & longitude<105 & latitude>=3.75 & latitude<=10 ~ 4,
    longitude>=105 & longitude<120 & latitude>=-15 & latitude<(-8.75) ~ 5,
    longitude>=105 & longitude<120 & latitude>=-8.75 & latitude<(-2.5) ~ 6,
    longitude>=105 & longitude<120 & latitude>=-2.5 & latitude<3.75 ~ 7,
    longitude>=105 & longitude<120 & latitude>=3.75 & latitude<=10 ~ 8,
    longitude>=120 & longitude<135 & latitude>=-15 & latitude<(-8.75) ~ 9,
    longitude>=120 & longitude<135 & latitude>=-8.75 & latitude<(-2.5) ~ 10,
    longitude>=120 & longitude<135 & latitude>=-2.5 & latitude<3.75 ~ 11,
    longitude>=120 & longitude<135 & latitude>=3.75 & latitude<=10 ~ 12,
    longitude>=135 & longitude<=150 & latitude>=-15 & latitude<(-8.75) ~ 13,
    longitude>=135 & longitude<=150 & latitude>=-8.75 & latitude<(-2.5) ~ 14,
    longitude>=135 & longitude<=150 & latitude>=-2.5 & latitude<3.75 ~ 15,
    longitude>=135 & longitude<=150 & latitude>=3.75 & latitude<=10 ~ 16,
  ))


df_10 <- air_temp_month_10 %>% 
  group_by(year,month,section) %>%  
  summarise(meanTemp = mean(air_temp))

air_temp_month_11 <- air_temp_month_11 %>% 
  mutate(section = case_when(
    longitude>=90 & longitude<105 & latitude>=-15 & latitude<(-8.75) ~ 1,
    longitude>=90 & longitude<105 & latitude>=-8.75 & latitude<(-2.5) ~ 2,
    longitude>=90 & longitude<105 & latitude>=-2.5 & latitude<3.75 ~ 3,
    longitude>=90 & longitude<105 & latitude>=3.75 & latitude<=10 ~ 4,
    longitude>=105 & longitude<120 & latitude>=-15 & latitude<(-8.75) ~ 5,
    longitude>=105 & longitude<120 & latitude>=-8.75 & latitude<(-2.5) ~ 6,
    longitude>=105 & longitude<120 & latitude>=-2.5 & latitude<3.75 ~ 7,
    longitude>=105 & longitude<120 & latitude>=3.75 & latitude<=10 ~ 8,
    longitude>=120 & longitude<135 & latitude>=-15 & latitude<(-8.75) ~ 9,
    longitude>=120 & longitude<135 & latitude>=-8.75 & latitude<(-2.5) ~ 10,
    longitude>=120 & longitude<135 & latitude>=-2.5 & latitude<3.75 ~ 11,
    longitude>=120 & longitude<135 & latitude>=3.75 & latitude<=10 ~ 12,
    longitude>=135 & longitude<=150 & latitude>=-15 & latitude<(-8.75) ~ 13,
    longitude>=135 & longitude<=150 & latitude>=-8.75 & latitude<(-2.5) ~ 14,
    longitude>=135 & longitude<=150 & latitude>=-2.5 & latitude<3.75 ~ 15,
    longitude>=135 & longitude<=150 & latitude>=3.75 & latitude<=10 ~ 16,
  ))

df_11 <- air_temp_month_11 %>% 
  group_by(year,month,section) %>%  
  summarise(meanTemp = mean(air_temp))

df <- rbind(df_10,df_11)


# select only the columns we wish to keep
df <- subset(df, select = c("year","month","section","meanTemp"))

# remove data that is duplicated. 
df <- df[!duplicated(df), ]

###################################################################
# add means to the "section"####
df_coll$year <- as.numeric(df_coll$year)
df_coll$month <- as.numeric(df_coll$month)
df_coll$section <- as.numeric(df_coll$section)

final_df <- df_coll %>% 
  mutate(mean_Air_Temp = case_when(
    (year == 2010 & month == 06 & section == 6) ~ 83.26881,
    (year == 2010 & month == 07 & section == 6) ~ 81.92102,
    (year == 2010 & month == 08 & section == 6) ~ 82.26683,
    (year == 2011 & month == 06 & section == 6) ~ 82.46824,
    (year == 2011 & month == 07 & section == 6) ~ 81.54429,
    (year == 2011 & month == 08 & section == 6) ~ 81.19775,
    
    (year == 2010 & month == 06 & section == 7) ~ 83.76042,
    (year == 2010 & month == 07 & section == 7) ~ 82.24087,
    (year == 2010 & month == 08 & section == 7) ~ 83.05614,
    (year == 2011 & month == 06 & section == 7) ~ 83.84140,
    (year == 2011 & month == 07 & section == 7) ~ 82.70522,
    (year == 2011 & month == 08 & section == 7) ~ 83.02953,
    
    (year == 2010 & month == 06 & section == 8) ~ 85.08822,
    (year == 2010 & month == 07 & section == 8) ~ 83.78226,
    (year == 2010 & month == 08 & section == 8) ~ 83.78701,
    (year == 2011 & month == 06 & section == 8) ~ 83.74804,
    (year == 2011 & month == 07 & section == 8) ~ 83.69563,
    (year == 2011 & month == 08 & section == 8) ~ 83.61275,
    
    (year == 2010 & month == 06 & section == 10) ~ 81.74605,
    (year == 2010 & month == 07 & section == 10) ~ 81.70000,
    (year == 2010 & month == 08 & section == 10) ~ 81.31765,
    (year == 2011 & month == 06 & section == 10) ~ 80.42105,
    (year == 2011 & month == 07 & section == 10) ~ 80.44911,
    (year == 2011 & month == 08 & section == 10) ~ 79.43655,
    
    (year == 2010 & month == 06 & section == 11) ~ 82.97083,
    (year == 2010 & month == 07 & section == 11) ~ 82.90743,
    (year == 2010 & month == 08 & section == 11) ~ 82.69141,
    (year == 2011 & month == 06 & section == 11) ~ 82.32911,
    (year == 2011 & month == 07 & section == 11) ~ 81.77465,
    (year == 2011 & month == 08 & section == 11) ~ 80.65878))

final_df <- subset(final_df, select = c("library_name","mean_Air_Temp")) %>% 
  dplyr::rename(sample = library_name)

#if you want to save the output--
write_csv(final_df,"./data/mean_air_temp_data/mean_air_temp_final.csv")  

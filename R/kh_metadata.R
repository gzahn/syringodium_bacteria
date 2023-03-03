# Kate's Metadata ####

##----------------
##SETUP
##----------------

library(tidyverse); packageVersion("tidyverse") # 1.3.2
library(geosphere); packageVersion("geosphere") # 1.15.18
library(rlist); packageVersion("rlist") # 0.4.6.2
library(dplyr); packageVersion("dplyr") # 1.0.10
library(purrr); packageVersion("purrr") # 0.3.4
library(stringr); packageVersion("stringr") # 1.4.1

##----------------
##CALL METADATA AND DATES
##----------------

# download csv files###
# collection dates by Wainwright (you need to make sure the location is right)
coll_dates <- read.csv("./metadata/syringodium_dates.csv")
metadata <- readRDS("./metadata/metadata_for_taxonomy.RDS")

##----------------
##COMMAND LINE: WGET URL
##----------------

# wget urls for each sea lvl station's csv in desired file system 
  # metadata directory for me

# navigate to desired directory to store 6 sealvl station csv 
  # uncomment & copy the code below, then paste into terminal:

# wget https://tidesandcurrents.noaa.gov/sltrends/data/550-003_meantrend.csv
# wget https://tidesandcurrents.noaa.gov/sltrends/data/550-005_meantrend.csv
# wget https://tidesandcurrents.noaa.gov/sltrends/data/550-007_meantrend.csv
# wget https://tidesandcurrents.noaa.gov/sltrends/data/550-009_meantrend.csv
# wget https://tidesandcurrents.noaa.gov/sltrends/data/550-014_meantrend.csv
# wget https://tidesandcurrents.noaa.gov/sltrends/data/550-017_meantrend.csv


# remove periods from the end of columns names manually
  #last two columns on confidence intervals have a period at the end of the title
  # (will find programmatic solution)


##----------------
##WRANGLE SEA LEVEL STATION
##----------------

# 6 Malay stations

## 550-003 Pulau Pinang, Malaysia
pulau_550_003 <- read_csv("./data/550-003_meantrend.csv")
pulau_550_003$Latitutude <- c("5.421667")
pulau_550_003$Longitude <- c("100.346667")
pulau_550_003$Name <- c("pulau")
pulau_550_003$sealvl_statid <- c("550_003")

## 550-005 Lumut, Malaysia
lumut_550_005 <- read_csv("./data/550-005_meantrend.csv")
lumut_550_005$Latitutude <- c("4.24")
lumut_550_005$Longitude <- c("100.613333")
lumut_550_005$Name <- c("lumut")
lumut_550_005$sealvl_statid <- c("550_005")

## 550-007 Pelabuhan Kelang, Malaysia
pelabuhan_550_007 <- read_csv("./data/550-007_meantrend.csv")
pelabuhan_550_007$Latitutude <- c("3.05")
pelabuhan_550_007$Longitude <- c("101.358333")
pelabuhan_550_007$Name <- c("pelabuhan")
pelabuhan_550_007$sealvl_statid <- c("550_007")

## 550-009 Tanjung Keling, Malaysia
tan_kel_550_009 <- read_csv("./data/550-009_meantrend.csv")
tan_kel_550_009$Latitutude <- c("2.215")
tan_kel_550_009$Longitude <- c("102.153333")
tan_kel_550_009$Name <- c("tan_kel")
tan_kel_550_009$sealvl_statid <- c("550_009")

## 550-014 Tanjung Gelang, Malaysia
tan_gel_550_014 <- read_csv("./data/550-014_meantrend.csv")
tan_gel_550_014$Latitutude <- c("3.975")
tan_gel_550_014$Longitude <- c("103.43")
tan_gel_550_014$Name <- c("tan_gel")
tan_gel_550_014$sealvl_statid <- c("550_014")

## 550-017 Cendering, Malaysia
cendering_550_017 <- read_csv("./data/550-017_meantrend.csv")
cendering_550_017$Latitutude <- c("5.265")
cendering_550_017$Longitude <- c("103.186667")
cendering_550_017$Name <- c("cendering")
cendering_550_017$sealvl_statid <- c("550_017")

# Create small metadata table with info about sea level stations
  ## name, id, id2, lat, long
sealvl_station <- data.frame(name=c("pulau", "lumut", "pelabuhan", "tan_kel", "tan_gel", "cendering"),
                             id=c("550-003", "550-005", "550-007", "550-009", "550-014", "550-017"),
                             id2=c("1594", "1595", "1591", "1593", "1589", "1592"),
                             lat=c("5.421667", "4.24", "3.05", "2.215", "3.975", "5.265"),
                             long=c("100.346667", "100.613333", "101.358333", "102.153333", "103.43", "103.186667"))

##----------------
##LINK SAMPLE & SEA LEVEL STATION
##----------------

# loop through samples & stations, closest station ~ lat + lon
  # prim_tabl is sample metadata & supp_tabl is weather station data
fun_closest_station <- function(prim_tabl, supp_tabl) {
  # Initialize new column for sea level
  prim_tabl$sealvl_statid <- NA
  # loop through metadata
  for (row in 1:nrow(prim_tabl)) {
    # save sample/lon/lat for current row
    samp_name <- prim_tabl[row, "sample"]
    samp_lon <- as.numeric(prim_tabl[row, "lon"])
    samp_lat <- as.numeric(prim_tabl[row, "lat"])
    # assign current lon/lat
    samp_loc <- c(samp_lon, samp_lat)
    
    # initializing closest distance from current sample
    curr_min <- 100000000 # very far away
    min_name <- ""
    
    # loop through 6 stations
    for (row2 in 1:nrow(supp_tabl)) {
      # save lon/lat for each site
      stat_lon <- as.numeric(supp_tabl[row2, "long"])
      stat_lat <- as.numeric(supp_tabl[row2, "lat"])
      # assign current lon/lat
      stat_loc <- c(stat_lon, stat_lat)

      # distHaversine() to get distance
      name <- supp_tabl[row2,2]
      dist <- distHaversine(samp_loc, stat_loc)
      if (dist < curr_min) {
        curr_min <- dist
        min_name <- name
      }
    }
    # add shortest-distance stationid to prim_tabl (metadata) sealvl_statid
    prim_tabl$sealvl_statid[row] <- min_name # fix it
  } 
  # return metadata table
  return(prim_tabl)
}

# call function to initialize object with sea level station id
metadata_site <- fun_closest_station(metadata, sealvl_station)

##----------------
##LINK SAMPLE & DATE
##----------------

# select rows that have a date for collection and no NA values
sam_dates <- coll_dates %>% 
  select(Library.Name, Collection_Date) %>% 
  filter(!is.na(Collection_Date) & !Library.Name == "SI")

# sort alphabetically so order matches metadata_site
sam_dates <- sam_dates[order(sam_dates$Library.Name), ]

# populate metadata_site$collection_date with sam_dates$Collection_dates by sample
  ## add a date column to metadata_site
metadata_site$collection_date <- NA

# loop through the metadata samples
for (row in 1:nrow(metadata_site)) {
  sample_md <- metadata_site[row, "sample"]
  # loop through the collection dates
  for (row2 in 1:nrow(sam_dates)) {
    sample_sd <- sam_dates[row2, "Library.Name"]
    # check to see if the sample names match
    if (sample_sd == sample_md) {
      # if the sample names match, add collection date to metadata
      metadata_site$collection_date[row] <- sam_dates$Collection_Date[row2]
    }    
  }
}

##----------------
##LINK SAMPLE & SEALVL ~ DATE
##----------------

# initiate column for monthly mean sea level (MMSL)
metadata_site$mmsl <- NA

# pipe select !na and set as var
metadata_site <- metadata_site %>% 
  filter(!is.na(collection_date))


# loop through samples to locate station id
for (row in 1:nrow(metadata_site)) {
  # select station id & sample
  stat_id <- metadata_site[row, "sealvl_statid"]
  if (stat_id == "550-003") {
    # fix and assign sample collection date by year and month: e.g. 2011_06
    dat <- str_split(metadata_site$collection_date[row], "-")
    dyear <- dat %>% map(1)
    # months are in "01" format, not '1' like station data
    dmonth <- dat %>% map(2)

    # loop through closest station to retrieve sea level (m)
    for (row2 in 1:nrow(pulau_550_003)) {
      smonth <- pulau_550_003$Month[row2]
      # if it's a single-number month, add a 0 to beginning
      if (nchar(smonth) < 2) {
        smonth <- paste0("0", smonth)
      }
      if (!is.na(dyear) && dyear == pulau_550_003$Year[row2] 
          & !is.na(dmonth) && dmonth == smonth) {
        # add sea level data to metadata_site
        metadata_site$mmsl[row] <- pulau_550_003$Monthly_MSL[row2]
      }
    }
  } else if (stat_id == "550-005") {
    # fix and assign sample collection date by year and month: e.g. 2011_06
    dat <- str_split(metadata_site$collection_date[row], "-")
    dyear <- dat %>% map(1)
    # these months are in "01" format, not '1' like station data
    dmonth <- dat %>% map(2)

    # loop through closest station to retrieve sea level (m)
    for (row2 in 1:nrow(lumut_550_005)) {
      smonth <- lumut_550_005$Month[row2]
      # if it's a single-number month, add a 0 to beginning
      if (nchar(smonth) < 2) {
        smonth <- paste0("0", smonth)
      }
      if (!is.na(dyear) && dyear == lumut_550_005$Year[row2] 
          & !is.na(dmonth) && dmonth == smonth) {
        # add sea level data to metadata_site
        metadata_site$mmsl[row] <- lumut_550_005$Monthly_MSL[row2]
      }
    }
  } else if (stat_id == "550-007") {
    # fix and assign sample collection date by year and month: e.g. 2011_06
    dat <- str_split(metadata_site$collection_date[row], "-")
    dyear <- dat %>% map(1)
    # these months are in "01" format, not '1' like station data
    dmonth <- dat %>% map(2)

    # loop through closest station to retrieve sea level (m)
    for (row2 in 1:nrow(pelabuhan_550_007)) {
      smonth <- pelabuhan_550_007$Month[row2]
      # if it's a single-number month, add a 0 to beginning
      if (nchar(smonth) < 2) {
        smonth <- paste0("0", smonth)
      }
      if (!is.na(dyear) && dyear == pelabuhan_550_007$Year[row2] 
          & !is.na(dmonth) && dmonth == smonth) {
        # add sea level data to metadata_site
        metadata_site$mmsl[row] <- pelabuhan_550_007$Monthly_MSL[row2]
      }
    }
  } else if (stat_id == "550-009") {
    # fix and assign sample collection date by year and month: e.g. 2011_06
    dat <- str_split(metadata_site$collection_date[row], "-")
    dyear <- dat %>% map(1)
    # these months are in "01" format, not '1' like station data
    dmonth <- dat %>% map(2)

    # loop through closest station to retrieve sea level (m)
    for (row2 in 1:nrow(tan_kel_550_009)) {
      smonth <- tan_kel_550_009$Month[row2]
      # if it's a single-number month, add a 0 to beginning
      if (nchar(smonth) < 2) {
        smonth <- paste0("0", smonth)
      }
      if (!is.na(dyear) && dyear == tan_kel_550_009$Year[row2] 
          & !is.na(dmonth) && dmonth == smonth) {
        # add sea level data to metadata_site
        metadata_site$mmsl[row] <- tan_kel_550_009$Monthly_MSL[row2]
      }
    }
  } else if (stat_id == "550-014") {
    # fix and assign date by year and month: e.g. 2011_06
    dat <- str_split(metadata_site$collection_date[row], "-")
    dyear <- dat %>% map(1)
    # these months are in "01" format, not '1' like station data
    dmonth <- dat %>% map(2)

    for (row2 in 1:nrow(tan_gel_550_014)) {
      smonth <- tan_gel_550_014$Month[row2]
      # if it's a single-number month, add a 0 to beginning
      if (nchar(smonth) < 2) {
        smonth <- paste0("0", smonth)
      }
      if (!is.na(dyear) && dyear == tan_gel_550_014$Year[row2] 
          & !is.na(dmonth) && dmonth == smonth) {
        # add sea level data to metadata_site
        metadata_site$mmsl[row] <- tan_gel_550_014$Monthly_MSL[row2]
      }
    }
  } else if (stat_id == "550-017") {
    dat <- str_split(metadata_site$collection_date[row], "-")
    # dat <- !is.na(dat)
    dyear <- dat %>% map(1)
    # these months are in "01" format, not '1' like station data
    dmonth <- dat %>% map(2)

    for (row2 in 1:nrow(cendering_550_017)) {
      smonth <- cendering_550_017$Month[row2]
      # if it's a single-number month, add a 0 to beginning
      if (nchar(smonth) < 2) {
        smonth <- paste0("0", smonth)
      }
      if (!is.na(dyear) && dyear == cendering_550_017$Year[row2] 
          & !is.na(dmonth) && dmonth == smonth) {
        # add sea level data to metadata_site
        metadata_site$mmsl[row] <- cendering_550_017$Monthly_MSL[row2]
      }
    }
  }
}

##----------------
##PREPARE CSV
##----------------

# I will need to make some changes if we need 192 rows 
# bc I removed all NAs before the for loop
md_test <- metadata_site %>% 
  select(sample, mmsl)

# rename columns to make it easy for Geoff
colnames(md_test)[1] = "sample_id"
colnames(md_test)[2] = "mean_month_slvl"

# write csv
md_test %>% 
  write_csv("./metadata/kh_metadata.csv")

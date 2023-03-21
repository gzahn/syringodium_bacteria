# CARBON MONOXIDE from Indonesia - NOAA
#The metadata of Carbon Monoxide samples was found by navigating to the NOAA website and then executing:
# wget https://gml.noaa.gov/aftp/data/trace_gases/co/flask/surface/txt/co_bkt_surface-flask_1_ccgg_event.txt
#This is found in "co_bkt_surface-flask_1_ccgg_event.txt"

#The Library list of the samples of microbiome collected is named "syringodium_dates.csv"


#The intent of this script is to create two columns of data. 
# The first column will come from the file "syringodium_dates.csv" 
#   and it will contain the Library name for each sample.
# The second column will come from the file "co_bkt_surface-flask_1_ccgg_event.txt"  and it will contain 
#   the sampled CO data for the matching date of the library sample from the "syringodium_dates.csv"  file. 
# Since the air CO samples were taken from one location in Indonesia, The
#   GPS coordinates of the two datafiles are not being used, but the dates will be solely used. 


###Note: The list of dates includes data from June, July, August of 2010 and June and August of 2011. 
### However, the metadata only provides sample data from the year 2010 for each of the months. 
### I tried searching CO2, Methane and other gases, and 
#no data was available for the year 2011 passed the month of March. 

#0. Load necessary Libraries: 
library(tidyverse)
library(skimr)
library(dplyr)
library(data.table)

#1. Getting the dates and library names from the collected samples file.  
#Navigate to desired directory containing the metadata.
#load in the "syringodium_dates.csv" containing the Library name for each sample. 
samData <- read.csv(file='./data/syringodium_dates.csv')
names(samData) #a glimpse at the columns

#selecting the dates and library sample columns and saving them into a new dataframe
dfDateSample <- samData %>% select('Collection_Date', 'Library.Name')
dfDateSample
#

#filling the empty slots with NA
dfDateSample[dfDateSample == ''] <- NA
dfDateSample

#checking the complete cases 
complete.cases(dfDateSample)

#if they are in complete (NA), then remove them
dfDateSampleCleaned <- dfDateSample[complete.cases(dfDateSample),]
dfDateSampleCleaned

head(dfDateSampleCleaned)

#renaming column names: 
colnames(dfDateSampleCleaned) <- c("Collection_Date", "Library_Name")
head(dfDateSampleCleaned)

#Final Dataframe containing the list of collection dates and library name
dfDateSampleCleaned



#2. Getting the dates and CO sample value from the metadata file
#load in the "co_bkt_surface-flask_1_ccgg_event.txt" containing the sampled CO data
coData<- read.table('./data/co_bkt_surface-flask_1_ccgg_event.txt', header = TRUE)
coData

names(coData)
#filter rows by year 2010 and 2011 & months 6, 7, 8
coDataSelectedDates<- coData %>% 
  filter (year %in% c("2010", "2011") & month %in% c("6","7","8"))

#select the columns of interest 
coDataSelectedCols <- coDataSelectedDates %>% select('datetime', 'value')
coDataSelectedCols

#renaming the rows of interests, for better data management
#   Jun 2010
coDataSelectedCols$datetime <- 
  #renaming the month of June
  gsub(x = coDataSelectedCols$datetime,  # Which column to change?
                pattern = "2010-06.*",  # change from
                replacement = "2010-06")  # change to   
#   Jul 2010
coDataSelectedCols$datetime <- 
   gsub(x = coDataSelectedCols$datetime, pattern = "2010-07.*", replacement = "2010-07")

#   Aug 2010
coDataSelectedCols$datetime <- 
   gsub(x = coDataSelectedCols$datetime, pattern = "2010-08.*",replacement = "2010-08")

#replacing datetime col title for date: 
colnames(coDataSelectedCols) <- c("Collection_Date", "CO_Value")


#MONTHLY AVERAGE
#calculating the average of data by month
#   June
meCOJun<- coDataSelectedCols %>% 
  filter (Collection_Date %in% c("2010-06"))
meCOJun<-mean (meCOJun$CO_Value)
meCOJun

#   July
meCOJul <- coDataSelectedCols %>% 
  filter (Collection_Date %in% c("2010-07"))
meCOJul<-mean (meCOJul$CO_Value)
meCOJul

#   Aug
meCOAug <- coDataSelectedCols %>% 
  filter (Collection_Date %in% c("2010-08"))
meCOAug<-mean (meCOAug$CO_Value)
meCOAug

#r the matching date of the library sample from the "syringodium_dates.csv"  file
Month_CO_Mean <- c(meCOJun, meCOJul, meCOAug)
Month_CO_Mean

Collection_Date <- c("2010-06", "2010-07", "2010-08") 

dfCODataDateMean <- data.frame(Collection_Date, Month_CO_Mean)
dfCODataDateMean

#3. Comparing the dataframe of the samples and metadata using dates
#   of both dataframes to  
#   assign a value of CO to each sample.

#joining the library name/dates dataframe (Step 1) & CO/dates dataframe (Step 2)
joinedMetadata <- full_join(dfCODataDateMean,dfDateSampleCleaned,by= 'Collection_Date')
View(joinedMetadata)

coMetadataFinal <- joinedMetadata %>% select('Library_Name', 'Month_CO_Mean')
coMetadataFinal
  
#4. Save and Create a new .csv file containing the Library name and CO samples. 
##put code in R folder 
write.csv(coMetadataFinal, "./data/co_metadata.csv")

#END
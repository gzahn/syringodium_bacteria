###########################################################
## Arctic Oscillation and Pacific West metadata addition ##
##Software versions:                                     ##
##                     tidyverse v 1.3.2                 ##
##                     curl v 5.0.0                      ##
###########################################################

# load necessary packages 
library(tidyverse); packageVersion("tidyverse")
library(curl); packageVersion("curl")

# read in the data 
dat<-read_csv("./data/syringodium_dates.csv")

fdate<-dat %>% drop_na(Collection_Date)

unique_collection_dates<-fdate$Collection_Date %>% unique()

###################################
####### Arctic Oscillation  #######
###################################

# the URL for the arctic oscillation data 
AOurl<-"https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ocean/index/heat_content_index.txt"

# how to get the URL into text 
ao_raw<-curl(AOurl)

AO_lines<-readLines(ao_raw)

# now use str_squish and then make a data frame only containing the things that we care about 
AO_dataframe<-str_squish(AO_lines) %>% data.frame()

ao_list1<-AO_dataframe[-1,] %>% data.frame()
ao_list<-ao_list1[-1,] %>% data.frame()
ao_list2<-ao_list[-1,] %>% data.frame()


behlee_ao<-separate(ao_list2, 
         col = '.',  
         into = c('year', 'month', "p1", "p2", "p3"), 
         sep=" ") 
# get the 2010 data and the 2011 data 
ba_ao_2010<-behlee_ao %>% filter(year==2010) 
ba_ao_2011<- behlee_ao %>% filter(year==2011)
# combine them 
ba_ao_2010_2011<-rbind(ba_ao_2010, ba_ao_2011)

# the date data is a little different in the arctic oscillation dataframe than it is in the collection data, so lets format it so it is easier down the line 
ba_ao_2010_2011$Collection_Date<-paste(ba_ao_2010_2011$year, ba_ao_2010_2011$month, sep = "-0")
  

#the data is collected at three locations, so lets take the mean of the data to simplify it
  # the data values are not numeric, so we have to fix that 
ba_ao_2010_2011$p1 = as.numeric(ba_ao_2010_2011$p1)
ba_ao_2010_2011$p2 = as.numeric(ba_ao_2010_2011$p2)
ba_ao_2010_2011$p3 = as.numeric(ba_ao_2010_2011$p3)

ba_ao_avg<-mutate(ba_ao_2010_2011, 
                  mean_ao=rowMeans(select(ba_ao_2010_2011, 
                                          p1:p3), na.rm = TRUE))

# now we can filter collection_date by the dates from our samples 
ba_ao_avg_cd<-filter(ba_ao_avg, Collection_Date %in% unique_collection_dates)
behlee_ao_cd<-ba_ao_avg_cd %>% select(Collection_Date, mean_ao) %>% data.frame()


#################################
####### West Pacific ############
#################################

# this is the url for the west pacific data 
WPurl<-"ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/wp_index.tim"
con<-curl(WPurl)

# this makes it so it something that sort of makes sense to a human 
wp_raw<-readLines(con)
# we have to squish this so that we can make it a dataframe and separate it later 
wp_squish<-str_squish(wp_raw) %>% data.frame()
    # there are rows  in here that are meaningless to the data so we need to remove them 
wp_squish1<-wp_squish[-1,] %>% data.frame()
wp_squish2<-wp_squish1[-1,] %>% data.frame()
wp_squish3<-wp_squish2[-1,] %>% data.frame() 
wp_squish4<-wp_squish3[-1,] %>% data.frame()
wp_squish5<-wp_squish4[-1,] %>% data.frame()
wp_squish6<-wp_squish5[-1,] %>% data.frame()
wp_squish7<-wp_squish6[-1,] %>% data.frame()
wp_squish8<-wp_squish7[-1,] %>% data.frame()

# made the dataframe with only the dates we want 
behlee_wp<-separate(wp_squish8, 
                    col = '.',  
                    into = c('year', 'month', "wp"), 
                    sep=" ") 

    # filter for the years we are interested in 
ba_wp_2010<-behlee_wp %>% filter(year==2010)
ba_wp_2011<- behlee_wp %>% filter(year==2011)
    # put them together so we can look for the specific dates that are in the collection data 
ba_wp_2010_2011<-rbind(ba_wp_2010, ba_wp_2011) 
    # the dates in the collection data are formatted a little differently than ours, so we can fix that here 
ba_wp_2010_2011$Collection_Date<-paste(ba_wp_2010_2011$year, ba_wp_2010_2011$month, sep = "-0")
    # now we can filter for the collection dates 
ba_wp_cd<-filter(ba_wp_2010_2011, Collection_Date %in% unique_collection_dates)
    # we have some extra columns that are not necisary, so we can get rid of everything except collection date and wp 
behlee_wp_cd<-ba_wp_cd %>% select(Collection_Date, wp) %>% data.frame

# now that we have both arctic oscillation and west pacific data with the same collection dates, we can put them together
behlee_meta<-full_join(behlee_ao_cd, behlee_wp_cd)



###########################################################
####### Adding the metadata to the sample data ############
###########################################################

# we can join the new data frame to the collection data to get the metadata associated with the sample names 
behlee_ao_wp_loc<-full_join(fdate, behlee_meta)

# to make it easier to use with our phyloseq object, we can select just the new metadata and the sample names 
wp_ao_id<-behlee_ao_wp_loc %>% select(`Library Name`, mean_ao, wp)
behlee_wpao_id<-wp_ao_id %>% dplyr::rename("sample" = "Library Name")

# wp is a character at this point, so we need to change it to numeric so we can do actual stuff with it
behlee_wpao_id$wp<-as.numeric(behlee_wpao_id$wp)

# wp and ao are only relevant to our data if we consider them in tandem. 
  # we need to find where both ao and wp are positive because that is when there are associations with temperature anomalies in south Asia 
    # to do this I need to do match pattern to do if the first one is positive, and the second one is also positive then return true
behlee_wpao_id_pos_wp<-behlee_wpao_id %>% 
  mutate(positive_wp=case_when((wp) > 0 ~ TRUE, 
                               TRUE ~ FALSE))


anomalies_and_othe_stuff<-behlee_wpao_id_pos_wp %>% 
  mutate(positive_ao=case_when(mean_ao > 0 ~ TRUE,
                               TRUE~FALSE), 
         temp_anomalies=case_when(positive_ao & positive_wp == TRUE ~ TRUE, 
                                  TRUE ~ FALSE))

temp_anomalies<-anomalies_and_othe_stuff %>% select(sample, temp_anomalies) %>% data.frame()

# lets write it as a csv for easy access later 
write_csv(temp_anomalies, "./data/BA_temp_anomalies_metadata.csv")


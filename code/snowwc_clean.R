# Dana Meadows snow survey data
# Data-cleaning script
# Data downloaded as csv from http://cdec.water.ca.gov/dynamicapp/wsSensorData
# Query parameters: station_id = dan, sensor # = 3, start_date = 1900-01-01, end_date = 2021-12-31
# View as csv, save as text file, open with LibreOffice Calc, save as csv. 

library(tidyverse)

# Read raw data and clean
snowwc_dan <- read_csv(here::here("data", "raw", "danameadows_snowsurvey.csv")) %>% 
  rename_with(tolower) %>% 
  rename(date = "date time") %>% 
  filter(str_detect(date, "04-01$")) %>% 
  mutate(snowwc_cm = (as.numeric(value) * 2.54),
         year = year(date),
         station_id = tolower(station_id),
         snowwc_pave = (snowwc_cm/mean(snowwc_cm)) * 100) %>%  # expressed as percent of average snowwc_cm as measured over entire time series
  select(station_id, date, year, snowwc_cm, snowwc_pave)

# Dana Meadows snow survey data
# Data-checking script for dataset created by snowwc_clean.R

library(readr)
library(dplyr)

# Create dataset
source(here::here("code", "snowwc_clean.R"))

# Check that there are no missing years
yr <- snowwc_dan %>% pull(year)
min <- min(yr)
max <- max(yr)
seq_yr <- seq(min, max, by = 1)
identical(yr, seq_yr)
rm(yr, min, max, seq_yr)

# Check that snowwc_pave are within acceptable range and without any NAs
range(snowwc_dan$snowwc_pave)
snowwc_dan %>% filter(is.na(snowwc_pave)) %>% count()


# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 4/1/2020
# Last Modified: 4/28/2020
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)

#-----------------------------------------------------------------
# read in biological data from trap
list.files('analysis/data/raw_data/WDFW')

bio_df = excel_sheets('analysis/data/raw_data/WDFW/AllBioData.xlsx') %>%
  as.list() %>%
  rlang::set_names() %>%
  map_df(.id = NULL,
         .f = function(yr) {
           data_df = read_excel('analysis/data/raw_data/WDFW/AllBioData.xlsx',
                                yr)
           return(data_df)
         }) %>%
  mutate_at(vars(TrapDate),
            list(ymd))

# put bounds around years
# min_yr = 2011
# max_yr = 2019
min_yr = min(bio_df$Year)
max_yr = max(bio_df$Year)


# pull out PIT tag numbers
tag_list = bio_df %>%
  split(list(.$Year)) %>%
  map(.f = function(x) {
    x %>%
      select(TagID)
  })

# save tags to upload to PTAGIS

# for(yr in names(tag_list)) {
#   write_delim(tag_list[[yr]],
#               path = paste0('analysis/data/raw_data/tag_lists/UC_Sthd_Tags_', yr, '.txt'),
#               delim = '\n',
#               col_names = F)
# }

# just write the latest year
write_delim(tag_list[[as.character(max_yr)]],
            path = paste0('analysis/data/raw_data/tag_lists/UC_Sthd_Tags_', max_yr, '.txt'),
            delim = '\n',
            col_names = F)

# save biological data for later
write_rds(bio_df,
          path = paste0('analysis/data/derived_data/Bio_Data_', min_yr, '_', max_yr, '.rds'))



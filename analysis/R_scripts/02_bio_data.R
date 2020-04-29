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

# tag_id = "3D9.1C2D8AE469"

# grab some extra data that was left out
extra_bio = read_excel("analysis/data/raw_data/WDFW/Steelhead_PRD_BY2011.xlsx",
                       "Tagging Data") %>%
  clean_names() %>%
  mutate(Year = 2011) %>%
  select(Year,
         TagID = pi_ttags,
         TrapDate = date,
         Sex = sex,
         Origin = origin,
         ForkLength = fork,
         Age = age) %>%
  filter(!is.na(TagID)) %>%
  bind_rows(read_excel("analysis/data/raw_data/WDFW/Steelhead_PRD_BY2012.xlsx",
                       "2011 STHD Tagging Data") %>%
              clean_names() %>%
              mutate(Year = 2012) %>%
              select(Year,
                     TagID = pi_ttags_1,
                     TrapDate = date,
                     Sex = sex,
                     Origin = origin,
                     ForkLength = fork,
                     Age = age) %>%
              bind_rows(read_excel("analysis/data/raw_data/WDFW/Steelhead_PRD_BY2012.xlsx",
                                   "UnknownSthdPRD2011Tagging") %>%
                          clean_names() %>%
                          mutate(Year = 2012) %>%
                          select(Year,
                                 TagID = tag_id,
                                 TrapDate = tag_date,
                                 Origin = rear_type,
                                 ForkLength = length)))

bio_df %>%
  filter(Year %in% c(2011, 2012)) %>%
  select(Year, TagID) %>%
  anti_join(extra_bio)


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



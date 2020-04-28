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
         })

# put bounds around years
# min_yr = 2011
# max_yr = 2019
min_yr = min(bio_df$Year)
max_yr = max(bio_df$Year)

# bio_df = min_yr:max_yr %>%
#   as.list %>%
#   rlang::set_names() %>%
#   map_df(.id = 'Year',
#          .f = function(yr) {
#
#            data_df = read_excel(paste0('analysis/data/raw_data/WDFW/UC_Steelhead_', yr, '.xlsx'),
#                            sheet = 'Bio Data')
#
#            if(c('PIT (Dorsal)') %in% names(data_df)) {
#              data_df %<>%
#                rename(`Primary Tag` = `PIT (Dorsal)`)
#            }
#
#            if(c('PIT (Unknown)') %in% names(data_df)) {
#              data_df %<>%
#                rename(`Secondary Tag` = `PIT (Unknown)`)
#            }
#
#            if(c('PIT(Unknown)') %in% names(data_df)) {
#              data_df %<>%
#                rename(`Secondary Tag` = `PIT(Unknown)`)
#            }
#
#            if(c('Primary PIT Tag') %in% names(data_df)) {
#              data_df %<>%
#                rename(`Primary Tag` = `Primary PIT Tag`)
#            }
#
#            if(c('Secondary PIT Tag') %in% names(data_df)) {
#              data_df %<>%
#                rename(`Secondary Tag` = `Secondary PIT Tag`)
#            }
#
#            if(c('Double Tag') %in% names(data_df)) {
#              data_df %<>%
#                rename(`Secondary Tag` = `Double Tag`)
#            }
#
#            data_df %<>%
#              rename(TagID = `Primary Tag`,
#                     TagOther = `Secondary Tag`) %>%
#              select(TagID,
#                     TagOther,
#                     TrapDate = SurveyDate,
#                     Origin = `Origin(final)`,
#                     Sex = `Sex(final)`,
#                     Age = FinalAge,
#                     ForkLength,
#                     POH,
#                     Weight) %>%
#              mutate(Origin = recode(Origin,
#                                     'h' = 'H'),
#                     Sex = recode(Sex,
#                                  'f' = 'F',
#                                  'm' = 'M')) %>%
#              distinct()
#
#            if(class(data_df$TrapDate) == 'character') {
#              data_df %<>%
#                mutate_at(vars(TrapDate),
#                          list(~ excel_numeric_to_date(as.numeric(.)))) %>%
#                mutate_at(vars(TrapDate),
#                          list(as.POSIXct))
#            }
#
#            # fix sex of one tag, an age 2 hatchery fish with "unknown" sex
#            if(yr == 2014) {
#              data_df %<>%
#                mutate(Sex = if_else(TagID == '3D9.1C2E046EA3',
#                                     'M',
#                                     Sex))
#            }
#
#            return(data_df)
#
#          })




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



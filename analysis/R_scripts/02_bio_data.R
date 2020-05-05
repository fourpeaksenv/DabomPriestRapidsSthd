# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 4/1/2020
# Last Modified: 5/5/2020
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(WriteXLS)

#-----------------------------------------------------------------
# read in biological data from trap
list.files('analysis/data/raw_data/WDFW')

bio_df = excel_sheets('analysis/data/raw_data/WDFW/PRD_BiologicalData_BY11-BY15.xlsx')[-1] %>%
  as.list() %>%
  rlang::set_names() %>%
  map_df(.id = 'BroodYear',
         .f = function(yr) {
    data_df = read_excel('analysis/data/raw_data/WDFW/PRD_BiologicalData_BY11-BY15.xlsx',
                         yr)
    # if("PIT (Dorsal)" %in% names(data_df)) {
    #   data_df %<>%
    #     rename(`PIT (Pelvic)` = `PIT (Dorsal)`)
    # }
    return(data_df)
  }) %>%
  bind_rows(read_excel('analysis/data/raw_data/WDFW/Steelhead_PRD_BY2016_QCI.xlsx',
                       "BioData") %>%
              mutate(BroodYear = "BY16")) %>%
  bind_rows(read_excel('analysis/data/raw_data/WDFW/Steelhead_PRD_BY2017_FlatFile.xlsx',
                       1) %>%
              mutate(BroodYear = "BY17")) %>%
  bind_rows(read_excel('analysis/data/raw_data/WDFW/BY18 BioData.xlsx',
                       1) %>%
              mutate(BroodYear = "BY18")) %>%
  bind_rows(read_excel('analysis/data/raw_data/WDFW/BY19 BioData.xlsx',
                       1) %>%
              mutate(BroodYear = "BY19")) %>%
  mutate(record_id = 1:n()) %>%
  mutate(Year = paste0("20", str_remove(BroodYear, "^BY")),
         Year = as.numeric(Year)) %>%
  gather(tag_loc, TagID, matches('^PIT')) %>%
  filter(!is.na(TagID)) %>%
  select(record_id,
         BroodYear,
         Year,
         tag_loc,
         TagID,
         Species = `Species(final)`,
         TrapDate = SurveyDate,
         Sex = `Sex(final)`,
         Origin = `Origin(final)`,
         ForkLength,
         Age = `Age (scales)`,
         FinalAge,
         AdClip = `Ad-clip`,
         CWT = `CWT (Sn)`) %>%
  mutate(Age = str_replace(Age, '^r', 'R')) %>%
  arrange(Year, record_id, tag_loc)

# save as Excel file
bio_df %>%
  split(list(.$Year)) %>%
  WriteXLS('analysis/data/raw_data/WDFW/PRA_Sthd_BioData.xlsx')

# put bounds around years
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



# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 4/1/2020
# Last Modified: 11/12/2020
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(openxlsx)

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

# any duplicated tags?
bio_df %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  arrange(Year, TagID, TrapDate, record_id) %>%
  tabyl(Year)

#-----------------------------------------------------------------
# add 2020 bio data
bio_2020 = read_csv('analysis/data/raw_data/WDFW/NBD-2019-189-PRD 1.csv') %>%
  mutate(BroodYear = "BY20") %>%
  mutate(Year = paste0("20", str_remove(BroodYear, "^BY")),
         Year = as.numeric(Year)) %>%
  mutate(Species = "ST") %>%
  pivot_longer(cols = matches('^PIT'),
               names_to = "tag_loc",
               values_to = "TagID") %>%
  filter(!is.na(TagID)) %>%
  mutate(record_id = seq(from = max(bio_df$record_id) + 1,
                         by = 1,
                         length.out = n())) %>%
  select(record_id,
         BroodYear,
         Year,
         tag_loc,
         TagID,
         Species,
         TrapDate = `Event Date`,
         Sex,
         Origin = `Final Origin`,
         ForkLength = Length,
         Age = `Scale Age`,
         AdClip = `Adipose Clip`,
         FinClip = `Fin Clip`,
         CWT_snout = `CWT - Snout`,
         CWT_body = `CWT - Body`) %>%
  mutate(TrapDate = mdy_hm(TrapDate)) %>%
  mutate(Age = str_replace(Age, '^r', 'R')) %>%
  mutate(AdClip = if_else(!is.na(AdClip),
                          "AD",
                          NA_character_)) %>%
  mutate(Sex = recode(Sex,
                      "Female" = "F",
                      "Male" = "M")) %>%
  mutate(across(c(starts_with("CWT")),
                ~ if_else(!is.na(.), T, F))) %>%
  mutate(CWT = if_else(CWT_snout,
                       "SN",
                       if_else(CWT_body,
                               "BD",
                               NA_character_))) %>%
  mutate(age_split = str_split(Age, "\\."),
         FinalAge = map_dbl(age_split,
                            .f = function(x) {
                              sum(as.numeric(x[1]),
                                  as.numeric(x[2]))
                            })) %>%
  select(all_of(names(bio_df)))

# any duplicated tags?
bio_2020 %>%
  filter(TagID %in% TagID[duplicated(TagID)])

#-----------------------------------------------------------------
# add to overall list
bio_df %<>%
  bind_rows(bio_2020)

#-----------------------------------------------------------------
# save as Excel file
#-----------------------------------------------------------------
bio_df %>%
  split(list(.$Year)) %>%
  write.xlsx(file = 'analysis/data/derived_data/PRA_Sthd_BioData.xlsx')

#-----------------------------------------------------------------
# for tag lists
#-----------------------------------------------------------------
# put bounds around years
min_yr = min(bio_df$Year)
max_yr = max(bio_df$Year)


# pull out PIT tag numbers
tag_list = bio_df %>%
  split(list(.$Year)) %>%
  map(.f = function(x) {
    x %>%
      select(TagID) %>%
      distinct()
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



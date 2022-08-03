# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 4/1/2020
# Last Modified: 1/18/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(writexl)
library(here)

#-----------------------------------------------------------------
# read in biological data from trap
list.files(here('analysis/data/raw_data/WDFW'))

bio_raw = excel_sheets(here('analysis/data/raw_data/WDFW/PRD_BiologicalData_BY11-BY15.xlsx'))[-1] %>%
  as.list() %>%
  rlang::set_names() %>%
  map_df(.id = 'BroodYear',
         .f = function(yr) {
    data_df = read_excel(here('analysis/data/raw_data/WDFW/PRD_BiologicalData_BY11-BY15.xlsx'),
                         yr)
    return(data_df)
  }) %>%
  bind_rows(read_excel(here('analysis/data/raw_data/WDFW/Steelhead_PRD_BY2016_QCI.xlsx'),
                       "BioData") %>%
              mutate(BroodYear = "BY16")) %>%
  bind_rows(read_excel(here('analysis/data/raw_data/WDFW/Steelhead_PRD_BY2017_FlatFile.xlsx'),
                       1) %>%
              mutate(BroodYear = "BY17")) %>%
  bind_rows(read_excel(here('analysis/data/raw_data/WDFW/BY18 BioData.xlsx'),
                       1) %>%
              mutate(BroodYear = "BY18")) %>%
  bind_rows(read_excel(here('analysis/data/raw_data/WDFW/BY19 BioData.xlsx'),
                       1) %>%
              mutate(BroodYear = "BY19")) %>%
  mutate(record_id = 1:n()) %>%
  mutate(Year = paste0("20", str_remove(BroodYear, "^BY")),
         Year = as.numeric(Year))

bio_raw %>%
  select(BroodYear, record_id, starts_with("PIT")) %>%
  pivot_longer(cols = starts_with("PIT"),
               names_to = "tag_loc",
               values_to = "tag_id") %>%
  filter(!is.na(tag_id)) %>%
  filter(tag_loc == "PIT (Unknown)") %>%
  select(record_id) %>%
  left_join(bio_raw %>%
              select(BroodYear, record_id, starts_with("PIT"))) %>%
  filter(!is.na(`PIT (Pelvic)`) |
           !is.na(`PIT (Dorsal)`)) %>%
  as.data.frame()

bio_df = bio_raw %>%
  pivot_longer(cols = starts_with("PIT"),
               names_to = "tag_loc",
               values_to = "tag_id") %>%
  filter(!is.na(tag_id)) %>%
  mutate(tag_loc = if_else(tag_loc %in% c("PIT (Pelvic)", "PIT (Dorsal)"),
                           "tag_code",
                           "tag_other")) %>%
  pivot_wider(names_from = "tag_loc",
              values_from = "tag_id") %>%
  mutate(tag_code = if_else(is.na(tag_code) & !is.na(tag_other),
                            tag_other,
                            tag_code),
         tag_other = if_else(tag_code == tag_other,
                             NA_character_,
                             tag_other)) %>%
  select(record_id,
         brood_year = BroodYear,
         year = Year,
         # tag_loc,
         tag_code,
         tag_other,
         species = `Species(final)`,
         trap_date = SurveyDate,
         sex = `Sex(final)`,
         origin = `Origin(final)`,
         fork_length = ForkLength,
         age = `Age (scales)`,
         final_age = FinalAge,
         ad_clip = `Ad-clip`,
         cwt = `CWT (Sn)`) %>%
  mutate(age = str_replace(age, '^r', 'R')) %>%
  arrange(year, record_id, tag_code)

# any duplicated tags?
bio_df %>%
  filter(tag_code %in% tag_code[duplicated(tag_code)]) %>%
  arrange(year, tag_code, trap_date, record_id) %>%
  tabyl(year)

# fish with more than one tag?
bio_df %>%
  filter(record_id %in% record_id[duplicated(record_id)]) %>%
  tabyl(year)

bio_df %>%
  filter(!is.na(tag_other))

#-----------------------------------------------------------------
# add 2020 bio data
bio_2020 = read_csv(here('analysis/data/raw_data/WDFW/NBD-2019-189-PRD 1.csv')) %>%
  janitor::clean_names() %>%
  mutate(brood_year = "BY20") %>%
  mutate(year = paste0("20", str_remove(brood_year, "^BY")),
         year = as.numeric(year)) %>%
  mutate(species = "ST") %>%
  mutate(record_id = seq(from = max(bio_df$record_id) + 1,
                         by = 1,
                         length.out = n())) %>%
  select(record_id,
         brood_year,
         year,
         tag_code = pit_tag,
         species,
         trap_date = event_date,
         sex,
         origin = final_origin,
         fork_length = length,
         age = scale_age,
         ad_clip = adipose_clip,
         fin_clip = fin_clip,
         cwt_snout,
         cwt_body) %>%
  mutate(trap_date = mdy_hm(trap_date)) %>%
  mutate(age = str_replace(age, '^r', 'R')) %>%
  mutate(ad_clip = if_else(!is.na(ad_clip),
                          "AD",
                          NA_character_)) %>%
  mutate(sex = recode(sex,
                      "Female" = "F",
                      "Male" = "M")) %>%
  mutate(across(c(starts_with("CWT")),
                ~ if_else(!is.na(.), T, F))) %>%
  mutate(cwt = if_else(cwt_snout,
                       "SN",
                       if_else(cwt_body,
                               "BD",
                               NA_character_))) %>%
  mutate(age_split = str_split(age, "\\."),
         final_age = map_dbl(age_split,
                            .f = function(x) {
                              sum(as.numeric(x[1]),
                                  as.numeric(x[2]))
                            })) %>%
  select(any_of(names(bio_df)))

# any duplicated tags?
bio_2020 %>%
  filter(tag_code %in% tag_code[duplicated(tag_code)])

#-----------------------------------------------------------------
# add to overall list
bio_df %<>%
  bind_rows(bio_2020)

#-----------------------------------------------------------------
# add 2021 bio data
# pull data from XML file from PTAGIS
library(XML)
library(methods)
library(xml2)

# use locally downloaded copy
xml_file = here('analysis/data/raw_data/WDFW/CME-2020-192-PRD.XML')
# query it from PTAGIS
# xml_file = "https://api.ptagis.org/files/mrr/CME-2020-192-PRD.XML"

data <- xml2::read_xml(xml_file)
doc <- XML::xmlParse(data)
col_nms <- XML::xmlToDataFrame(nodes = getNodeSet(doc, "//DetailProjectDefinedField")) %>%
  as_tibble()
df <- XML::xmlToDataFrame(nodes = getNodeSet(doc, "//MRREvent")) %>%
  as_tibble()
names(df)[str_detect(names(df), "PDV")] <- col_nms$Label[match(names(df)[str_detect(names(df), "PDV")], col_nms$PDVColumn)]

bio_2021 <- df %>%
  clean_names() %>%
  mutate(across(event_date,
                ymd_hms)) %>%
  mutate(brood_year = "BY21") %>%
  mutate(year = paste0("20", str_remove(brood_year, "^BY")),
         year = as.numeric(year)) %>%
  mutate(species = "ST") %>%
  mutate(record_id = seq(from = max(bio_df$record_id) + 1,
                         by = 1,
                         length.out = n())) %>%
  # add some additional information (age) from another file
  full_join(read_csv(here('analysis/data/raw_data/WDFW/CME-2020-192-PRD correct columns for Kevin final 2-10-21.csv')) %>%
              janitor::clean_names() %>%
              select(pit_tag,
                     tag_other = second_pit_tag,
                     final_origin,
                     scale_age) %>%
              distinct()) %>%
  mutate(across(length,
                as.numeric)) %>%
  # mutate(scale_age = NA_character_) %>%
  select(record_id,
         brood_year,
         year,
         tag_code = pit_tag,
         species,
         trap_date = event_date,
         sex,
         origin = final_origin,
         fork_length = length,
         age = scale_age,
         ad_clip = adipose_clip,
         fin_clip = fin_clip,
         cwt_snout = cwt_s,
         cwt_body,
         conditional_comments,
         srr = species_run_rear_type) %>%
  mutate(srr_origin = str_extract(srr, "[A-Z]")) %>%
  # mutate(trap_date = mdy_hm(trap_date)) %>%
  mutate(age = str_replace(age, '^r', 'R')) %>%
  mutate(ad_clip = if_else(!is.na(ad_clip),
                           "AD",
                           NA_character_)) %>%
  mutate(sex = recode(sex,
                      "Female" = "F",
                      "Male" = "M")) %>%
  mutate(across(c(starts_with("cwt")),
                ~ if_else(!is.na(.), T, F))) %>%
  mutate(cwt = if_else(cwt_snout,
                       "SN",
                       if_else(cwt_body,
                               "BD",
                               NA_character_))) %>%
  mutate(age_split = str_split(age, "\\."),
         final_age = map_dbl(age_split,
                             .f = function(x) {
                               sum(as.numeric(x[1]),
                                   as.numeric(x[2]))
                             })) %>%
  # select(any_of(names(bio_df)))
  select(any_of(names(bio_df)), everything())

# fix a couple bits of information after checking with Ben Truscott
bio_2021 %<>%
  mutate(origin = if_else(tag_code == "3DD.003DA287CE",
                          "H",
                          origin))

bio_2021 %<>%
  filter(!(tag_code == "3DD.003DA2859A" & str_detect(conditional_comments, "AI")))

bio_2021 %<>%
  filter(!(tag_code == "3DD.003DA28690" & srr == "32W"))


# any duplicated tags?
bio_2021 %>%
  filter(tag_code %in% tag_code[duplicated(tag_code)]) %>%
  as.data.frame()

# Origin doesn't match SRR origin
bio_2021 %>%
  filter(origin != srr_origin) %>%
  pull(tag_code) %>%
  unique() %>%
  paste(collapse = ", ")

bio_2021 %>%
  filter(is.na(final_age)) %>%
  tabyl(age_split, age)




names(bio_df)[!names(bio_df) %in% names(bio_2021)]
names(bio_2021)[!names(bio_2021) %in% names(bio_df)]

#-----------------------------------------------------------------
# add to overall list
bio_df %<>%
  bind_rows(bio_2021 %>%
              select(any_of(names(bio_df))))

#-----------------------------------------------------------------
# add 2022

bio_2022 <- read_csv(here("analysis/data/raw_data/WDFW",
                          "2021 OLAFT Steelhead CME-2021-183-PRD final.csv")) %>%
  janitor::clean_names() %>%
  mutate(brood_year = "BY22") %>%
  mutate(year = paste0("20", str_remove(brood_year, "^BY")),
         year = as.numeric(year)) %>%
  mutate(species = "ST") %>%
  mutate(record_id = seq(from = max(bio_df$record_id) + 1,
                         by = 1,
                         length.out = n())) %>%
  select(record_id,
         brood_year,
         year,
         tag_code = pit_tag,
         species,
         trap_date = event_date,
         sex,
         origin = final_origin,
         fork_length = length,
         age = scale_age,
         ad_clip = adipose_clip,
         fin_clip = fin_clip,
         cwt_snout = cwt_s,
         cwt_body) %>%
  mutate(trap_date = mdy_hm(trap_date)) %>%
  mutate(age = str_replace(age, '^r', 'R')) %>%
  mutate(ad_clip = if_else(!is.na(ad_clip),
                           "AD",
                           NA_character_)) %>%
  mutate(sex = recode(sex,
                      "Female" = "F",
                      "Male" = "M")) %>%
  mutate(across(c(starts_with("CWT")),
                ~ if_else(!is.na(.), T, F))) %>%
  mutate(cwt = if_else(cwt_snout,
                       "SN",
                       if_else(cwt_body,
                               "BD",
                               NA_character_))) %>%
  mutate(age_split = str_split(age, "\\."),
         final_age = map_dbl(age_split,
                             .f = function(x) {
                               sum(as.numeric(x[1]),
                                   as.numeric(x[2]))
                             })) %>%
  select(any_of(names(bio_df)))

#-----------------------------------------------------------------
# add to overall list
bio_df %<>%
  bind_rows(bio_2022 %>%
              select(any_of(names(bio_df))))


#-----------------------------------------------------------------
# save as Excel file
#-----------------------------------------------------------------
bio_df %>%
  split(list(.$year)) %>%
  write_xlsx(path = here('analysis/data/derived_data',
                         'PRA_Sthd_BioData.xlsx'))

#-----------------------------------------------------------------
# for tag lists
#-----------------------------------------------------------------
# put bounds around years
min_yr = min(bio_df$year)
max_yr = max(bio_df$year)


# pull out PIT tag numbers
tag_list = bio_df %>%
  split(list(.$year)) %>%
  map(.f = function(x) {
    x %>%
      pivot_longer(cols = starts_with("tag"),
                   names_to = "source",
                   values_to = "tag_code") %>%
      filter(!is.na(tag_code)) %>%
      select(tag_code)
  })

# save tags to upload to PTAGIS

# for(yr in names(tag_list)) {
#   write_delim(tag_list[[yr]],
#               file = here('analysis/data/raw_data/tag_lists',
#                           paste0('UC_Sthd_Tags_', yr, '.txt')),
#               delim = '\n',
#               col_names = F)
# }

# just write the latest year
write_delim(tag_list[[as.character(max_yr)]],
            file = here('analysis/data/raw_data/tag_lists',
                        paste0('UC_Sthd_Tags_', max_yr, '.txt')),
            delim = '\n',
            col_names = F)

# save biological data for later
write_rds(bio_df,
          file = here('analysis/data/derived_data',
                      paste0('Bio_Data_', min_yr, '_', max_yr, '.rds')))



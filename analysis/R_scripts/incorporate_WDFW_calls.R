# Author: Kevin See
# Purpose: move old WDFW PITcleanr calls into new format
# Created: 12/2/21
# Last Modified: 12/2/21
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(lubridate)
library(janitor)
library(readxl)
library(magrittr)
library(here)

#-----------------------------------------------------------------
# set year
yr = 2015


# load PITcleanr prepped data
load(here('analysis/data/derived_data/PITcleanr',
          paste0('UC_Steelhead_', yr, '.rda')))

wdfw_obs = read_excel(here("analysis/data/derived_data/WDFW",
                           paste0("UC_Steelhead_", yr, ".xlsx"))) %>%
  clean_names() %>%
  mutate(across(ends_with("proc_status"),
                as.logical)) %>%
  rename(tag_code = tag_id,
         start_date = trap_date,
         min_det = obs_date,
         max_det = last_obs_date) %>%
  mutate(across(c(min_det, max_det),
                ymd_hms),
         across(start_date,
                ymd)) %>%
  select(tag_code, min_det, max_det,
         node, node_order,
         user_proc_status) %>%
  # early years coded LMRB0 as LMR (single array)
  mutate(across(node,
                str_replace,
                "LMR$",
                "LMRB0"),
         node_order = if_else(node == "LMRB0",
                              min(node_order[node == "LMRB0"]),
                              node_order),
         across(node,
                str_replace,
                "BelowJD1",
                "JDA"))

tag_paths = wdfw_obs %>%
  filter(user_proc_status) %>%
  select(tag_code, node, node_order) %>%
  distinct() %>%
  group_by(tag_code) %>%
  nest() %>%
  mutate(wdfw_keep_path = map_chr(data,
                                  .f = function(x) {
                                    x %>%
                                      arrange(node_order) %>%
                                      pull(node) %>%
                                      paste(collapse = " ")
                                  })) %>%
  select(-data) %>%
  ungroup() %>%
  full_join(prepped_ch %>%
              filter(auto_keep_obs) %>%
              select(tag_code, node, node_order) %>%
              distinct() %>%
              arrange(tag_code, node_order) %>%
              group_by(tag_code) %>%
              nest() %>%
              mutate(pitcl_keep_path = map_chr(data,
                                               .f = function(x) {
                                                 paste(x$node, collapse = " ")
                                               })) %>%
              select(-data) %>%
              ungroup())

tag_paths %>%
  filter(wdfw_keep_path != pitcl_keep_path)


tags_to_fix = prepped_ch %>%
  filter(is.na(user_keep_obs)) %>%
  left_join(wdfw_obs %>%
              select(tag_code,
                     node,
                     user_proc_status) %>%
              distinct()) %>%
  left_join(tag_paths) %>%
  filter(wdfw_keep_path != pitcl_keep_path) %>%
  select(tag_code, ends_with("_path")) %>%
  distinct()

# some tags have WEA in older file, and WEH in new file
wells_tags = tags_to_fix %>%
  filter((str_detect(wdfw_keep_path, "WEA$") &
            str_detect(pitcl_keep_path, "WEH")) |
           (str_detect(wdfw_keep_path, "WVT") &
              str_detect(pitcl_keep_path, "WEA$")))
tags_to_fix %<>%
  anti_join(wells_tags)

new_prepped = prepped_ch %>%
  anti_join(tags_to_fix) %>%
  mutate(user_keep_obs = auto_keep_obs) %>%
  bind_rows(prepped_ch %>%
              inner_join(tags_to_fix) %>%
              rowwise() %>%
              mutate(user_keep_obs = if_else(str_detect(wdfw_keep_path, node),
                                             T, F)) %>%
              ungroup()) %>%
  select(-ends_with("keep_path"))


tag_paths_2 = tag_paths %>%
  select(tag_code,
         wdfw_keep_path) %>%
  full_join(new_prepped %>%
              filter(user_keep_obs) %>%
              select(tag_code, node, node_order) %>%
              distinct() %>%
              arrange(tag_code, node_order) %>%
              group_by(tag_code) %>%
              nest() %>%
              mutate(pitcl_keep_path = map_chr(data,
                                               .f = function(x) {
                                                 paste(x$node, collapse = " ")
                                               })) %>%
              select(-data) %>%
              ungroup())

tag_paths_2 %>%
  filter(!tag_code %in% wells_tags$tag_code) %>%
  filter(wdfw_keep_path != pitcl_keep_path) %>%
  filter(str_detect(pitcl_keep_path, "WEH", negate = T)) %>%
  filter(!(str_detect(wdfw_keep_path, "WVT") & str_detect(pitcl_keep_path, "WEA$"))) %>%
  filter(str_detect(wdfw_keep_path, "ENS", negate = T),
         str_detect(wdfw_keep_path, "ENM", negate = T)) %>%
  filter(str_detect(pitcl_keep_path, "LWEA0", negate = T)) %>%
  filter(str_detect(pitcl_keep_path, "TWISPW", negate = T))

# overwrite saved files
prepped_ch = new_prepped
save(parent_child, configuration, start_date, bio_df, prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('UC_Steelhead_', yr, '.rda')))


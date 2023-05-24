# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 4/27/20
# Last Modified: 11/7/22
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(writexl)
library(here)

library(PITcleanr)

# library(XML)
# library(methods)
# library(xml2)
#
# #-----------------------------------------------------------------
# # function to read in XML file and transform into tibble
# readXMLTagFile = function(file_nm,
#                           file_path = NULL) {
#
#   if(!is.null(file_path)) {
#     data <- xml2::read_xml(paste(file_path, file_nm, sep = "/"))
#   } else {
#     data <- xml2::read_xml(paste("https://api.ptagis.org/files/mrr",
#                                  file_nm,
#                                  sep = "/"))
#   }
#   doc <- XML::xmlParse(data)
#   col_nms <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//DetailProjectDefinedField")) %>%
#     as_tibble()
#   df <- XML::xmlToDataFrame(nodes = XML::getNodeSet(doc, "//MRREvent")) %>%
#     as_tibble()
#   names(df)[str_detect(names(df), "PDV")] <- col_nms$Label[match(names(df)[str_detect(names(df), "PDV")], col_nms$PDVColumn)]
#
#   df <- df |>
#     janitor::clean_names() |>
#     dplyr::mutate(
#       dplyr::across(
#         sequence_number,
#         as.numeric
#       )
#     )
#
#   if(sum(str_detect(names(df), "date")) > 0) {
#     df <- df |>
#       mutate(across(contains("date"),
#                     lubridate::ymd_hms))
#   }
#
#
#   return(df)
# }
#
# readTagFile = function(file_nm,
#                        file_path = NULL,
#                        return_type = c("tag",
#                                        "meta",
#                                        "both")) {
#
#   return_type = match.arg(return_type)
#
#   if(!is.null(file_path)) {
#     my_content = paste(file_path,
#                        file_nm,
#                        sep = "/")
#   } else {
#     my_content <- httr::GET(paste("https://api.ptagis.org/files/mrr",
#                                   file_nm,
#                                   sep = "/")) |>
#       httr::warn_for_status() |>
#       httr::content("text",
#                     encoding = "UTF-8")
#   }
#
#   first_tag_row <-
#     my_content |>
#     readr::read_tsv(col_names = F,
#                     trim_ws = T,
#                     show_col_types = FALSE) |>
#     summarize(row_num = which(str_detect(X1, "^1 "))) |>
#     pull(row_num)
#
#   tag_data <- my_content |>
#     readr::read_tsv(col_names = F,
#                     skip = first_tag_row - 1,
#                     trim_ws = T,
#                     show_col_types = FALSE) |>
#     mutate(split_text = str_split(X1, "  ")) |>
#     mutate(id = map_chr(split_text,
#                         extract,
#                         1),
#            mutate(across(id,
#                          ~ as.numeric(.) |>
#                            suppressWarnings())),
#            pit_tag = map_chr(split_text,
#                              extract,
#                              2),
#            comments = map_chr(split_text,
#                               extract,
#                               12)
#     )
#
#   dates <- tag_data |>
#     dplyr::filter(is.na(id),
#                   stringr::str_detect(X1, "^CLOSE DATE",
#                                       negate = T)) |>
#     dplyr::select(X1) |>
#     dplyr::mutate(grp_num = stringr::str_split(X1, "=", simplify = T)[,1],
#                   event_date = stringr::str_split(X1, "=", simplify = T)[,2]) |>
#     dplyr::mutate(
#       dplyr::across(grp_num,
#                     ~ stringr::str_remove(., "^V")),
#       dplyr::across(event_date,
#                     ~ lubridate::mdy_hm(.))
#     ) |>
#     dplyr::select(grp_num,
#                   event_date) |>
#     dplyr::mutate(
#       dplyr::across(event_date,
#                     ~ lubridate::floor_date(., unit = "days")))
#
#   tag_data <- tag_data |>
#     filter(!is.na(id),
#            !is.na(pit_tag)) |>
#     mutate(comments_split = str_split(comments, "\\|"),
#            srr = map_chr(comments_split,
#                          extract,
#                          1),
#            grp_num = stringr::str_sub(srr, -2),
#            across(srr,
#                   ~ str_sub(., 1, 3)),
#            conditional_comments = map_chr(comments_split,
#                                           extract,
#                                           2),
#            text_comments = map_chr(comments_split,
#                                    extract,
#                                    3)) |>
#     mutate(event_type = if_else(str_detect(conditional_comments, "RE"),
#                                 "Recapture",
#                                 "Mark")) |>
#     left_join(dates,
#               by = join_by(grp_num)) |>
#     select(sequence_number = id,
#            pit_tag,
#            species_run_rear_type = srr,
#            event_date,
#            event_type,
#            conditional_comments,
#            text_comments)
#
#
#   meta_data <- my_content |>
#     readr::read_delim(delim = ":",
#                       trim_ws = T,
#                       n_max = first_tag_row - 1,
#                       col_names = c("name",
#                                     "value"),
#                       show_col_types = FALSE) |>
#     dplyr::filter(str_detect(name, "- - -", negate = T)) |>
#     tidyr::pivot_wider() |>
#     janitor::clean_names() |>
#     mutate(across(contains("date"),
#                   lubridate::mdy_hm))
#
#   if(return_type == "tag") {
#     return(tag_data)
#   } else if(return_type == "meta") {
#     return(meta_data)
#   } else {
#     return(list("tag" = tag_data,
#                 "meta" = meta_data))
#   }
# }

#-----------------------------------------------------------------
# where do tagging files live?
# an example of locally stored tagging files
data_path = paste0("T:/DFW-Team FP Upper Columbia Escapement - General/",
                   "UC_Sthd/inputs/PTAGIS/",
                   "Priest_Tagging_Files")

list.files(data_path)
xml_files = list.files(data_path)[str_detect(list.files(data_path), "xml$")]

all_tag_df <- tibble(file_nm = xml_files) |>
  mutate(tag_file = map(file_nm,
                        .f = function(x) {
                          queryMRRDataFile(
                            file_nm = x,
                            file_path = data_path)
                          })) |>
  unnest(tag_file)

names(all_tag_df)
nrow(all_tag_df)
tabyl(all_tag_df, file_nm) |>
  adorn_totals()


#-----------------------------------------------------------------
# now query all files from PTAGIS
ptagis_file_nms <- read_excel(paste(data_path,
                                    "PRD PTAGIS File Names and Information.xlsx",
                                    sep = "/"),
                              1) |>
  pull(1)

# tmp <- queryMRRDataFile(file_nm = ptagis_file_nms[9])
# names(tmp)
#
# tmp
# tabyl(tmp, species_run_rear_type)

all_tag_df <- tibble(file_nm = ptagis_file_nms) |>
  mutate(file_format = if_else(str_detect(file_nm,
                                        ".xml$"),
                             "xml",
                             "txt")) |>
  mutate(tag_file = map(file_nm,
                        .f = function(x) {
                          queryMRRDataFile(x)
                        })) |>
  unnest(tag_file)

tabyl(all_tag_df,
      species_run_rear_type) |>
  adorn_totals()

all_tag_df |>
  filter(str_detect(species_run_rear_type,
                    "[:alpha:]",
                    negate = T)) |>
  select(file_nm) |>
  distinct()


#-----------------------------------------------------------------
# compare with older tag lists
tag_list_path = "T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/inputs/PTAGIS/tag_lists"
tag_list_org <- tibble(tag_files = list.files(tag_list_path)) |>
  mutate(spawn_year = str_extract(tag_files,
                                  "[:digit:]+"),
         tag_df = map(tag_files,
                      .f = function(x) {
                        read_delim(paste(tag_list_path,
                                         x,
                                         sep = "/"),
                                   delim = "\n",
                                   col_names = "pit_tag")
                      })) |>
  unnest(tag_df) |>
  select(spawn_year, pit_tag)


all_tag_df |>
  select(-spawn_year) |>
  filter(str_detect(species_run_rear_type, "^3")) |>
  anti_join(tag_list_org) |>
  select(pit_tag:event_date, event_type,
         contains("comments"))
  tabyl(species_run_rear_type)

tag_list_org |>
  anti_join(all_tag_df |>
              select(pit_tag)) |>
  tabyl(spawn_year)

# try to pull out possible second tags
all_tag_df |>
  filter(str_detect(text_comments,
                    "SAME AS") |
           str_detect(text_comments,
                      "Same as") |
           str_detect(text_comments,
                      "3[:alnum:][:alnum:]\\.")) |>
  select(pit_tag, event_date, text_comments) |>
  mutate(across(text_comments,
                ~ str_remove(., "\\#"))) |>
  mutate(second_tag = str_extract(text_comments, "[:space:]3[:alnum:][:alnum:]\\.[:alnum:]+"),
         second_tag = if_else(is.na(second_tag),
                              str_extract(text_comments, "3D[:alnum:]\\. [:alnum:]+"),
                              second_tag),
         second_tag = if_else(is.na(second_tag),
                              paste0("3D9.", str_remove(text_comments, "SAME AS")),
                              second_tag)) |>
  mutate(across(second_tag,
                ~ str_remove(., "[:space:]+")),
         across(second_tag,
                ~ str_to_upper(.))) |>
  mutate(in_org_tags = if_else(second_tag %in% tag_list_org$pit_tag,
                               T, F)) |>
  # as.data.frame()
  select(pit_tag = second_tag) -> dbl_tags

tag_list_org |>
  anti_join(all_tag_df |>
              select(pit_tag) |>
              bind_rows(dbl_tags)) |>
  tabyl(spawn_year)


#-----------------------------------------------------------------
# decode conditional comments
cond_comm_codes <- read_csv(
  paste0("T:/DFW-Team FP Upper Columbia Escapement - General/",
         "UC_Sthd/inputs/PTAGIS/",
         "Glossary_ConditionalComment_ValidationCodes.csv")) |>
  clean_names()


all_tag_df |>
  select(#pit_tag,
    conditional_comments) |>
  # sample_n(1000) |>
  distinct() |>
  mutate(tmp = map(conditional_comments,
                   .f = function(x) {
                     str_split(x, "[:space:]") |>
                       unlist() |>
                       as_tibble() |>
                       rename(code = value)
                   })) |>
  select(-conditional_comments) |>
  unnest(tmp) |>
  # tabyl(code)
  left_join(cond_comm_codes |>
              select(-definition)) |>
  select(code, name) |> distinct()

all_tag_df |>
  filter(str_detect(conditional_comments, "DB")) |>
  as.data.frame() |> head()



tabyl(all_tag_df,
      species_run_rear_type,
      # migratory_yr)
      migration_year)
tabyl(all_tag_df,
      event_date,
      type) |>
  head()

all_tag_df |>
  group_by(file_nm,
           migration_year,
           type) |>
  summarize(across(event_date,
                   n_distinct),
            .groups = "drop") |>
  filter(event_date > 1)

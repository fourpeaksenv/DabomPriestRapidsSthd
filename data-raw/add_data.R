# Author: Kevin See
# Purpose: read in age table
# Created: 11/30/22
# Last Modified: 11/30/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(usethis)
library(janitor)

#-----------------------------------------------------------------
age_table = read_excel(here("inst/extdata",
                            "Age_Cheat_Sheet.xlsx"),
                       1,
                       skip = 1) |>
  clean_names() %>%
  rename(run = x3,
         age = x4,
         total_age = x9,
         notes = x10) %>%
  mutate(across(winter_spent_in_gravel_prior_to_hatching:total_age,
                as.numeric))

# save for use within "package"
use_data(age_table,
         overwrite = T)

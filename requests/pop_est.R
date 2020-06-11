# Author: Kevin See
# Purpose: Compile population level estimates for all years
# Created: 6/11/20
# Last Modified: 6/11/20
# Notes: Mark Miller requested population scale estimates for all years

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)

#-----------------------------------------------------------------
pop_est = 2011:2019 %>%
  as.list %>%
  rlang::set_names() %>%
  map_df(.id = 'SpawnYear',
         .f = function(yr) {
    read_excel(paste0('outgoing/estimates/PRA_Steelhead_', yr, '_20200604.xlsx'),
               "Population Escapement")
  }) %>%
  mutate_at(vars(SpawnYear),
            list(as.numeric)) %>%
  filter(Population %in% c('Wenatchee',
                           'Entiat',
                           'Methow',
                           'Okanogan')) %>%
  mutate(Origin = factor(Origin,
                         levels = c('Natural', 'Hatchery')),
         Population = factor(Population,
                             levels = c('Wenatchee',
                                        'Entiat',
                                        'Methow',
                                        'Okanogan'))) %>%
  arrange(SpawnYear, Population, Origin)

pop_est %>%
  write_csv("outgoing/other/UC_Sthd_Pop_Est.csv")

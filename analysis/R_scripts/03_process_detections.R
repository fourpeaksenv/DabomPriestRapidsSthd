# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 4/1/20
# Last Modified: 4/1/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(readxl)
library(magrittr)

#-----------------------------------------------------------------
# load configuration and site_df data
load('analysis/data/derived_data/site_config.rda')

# which spawn year are we dealing with?
yr = 2019

# start date is May 1
start_date = paste0(yr, '0501')

# build parent-child table
parent_child = createParentChildDf(site_df,
                                   configuration,
                                   startDate = start_date)

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
observations = read_csv(paste0('analysis/data/raw_data/PTAGIS/Tumwater_Chinook_', yr, '.csv'))

if(yr == 2019) {
  # 4 tags were collected at Tumwater and released into the White River above detection site, so should not be included in model
  observations %<>%
    filter(! `Tag Code` %in% c('3DD.00777D1003',
                               '3DD.00777D364F',
                               '3DD.0077B5E7BF',
                               '3DD.0077D002A6'))
}

# deal with some double tagged fish
bio_df = read_rds('analysis/data/derived_data/Bio_Data_2008_2019.rds') %>%
  filter(Year == yr)

dbl_tag = bio_df %>%
  filter(!is.na(TagOther))

if(nrow(dbl_tag) > 0) {
  observations %<>%
    left_join(dbl_tag %>%
                select(`Tag Code` = TagOther,
                       TagID)) %>%
    # filter(!is.na(TagID))
    mutate(`Tag Code` = if_else(!is.na(TagID),
                                TagID,
                                `Tag Code`)) %>%
    select(-TagID)
}

# # remove detections at UWE
# observations %<>%
#   filter(`Event Site Code Value` != 'UWE')

# process those observations with PITcleanr, using Tumwater-specific function
proc_list = processCapHist_TUM(start_date = start_date,
                               configuration = configuration,
                               parent_child = parent_child,
                               observations = observations,
                               last_obs_date = paste0(yr, "0930"),
                               truncate = T,
                               save_file = T,
                               file_name = paste0('outgoing/PITcleanr/TUM_Chinook_', yr, '.xlsx'))

# save some stuff
save(yr, start_date, parent_child, proc_list,
     file = paste0('analysis/data/derived_data/PITcleanr/TUM_Chinook_', yr, '.rda'))

#-----------------------------------------------------------------
# Read in the file returned by the Yakima
#-----------------------------------------------------------------
load(paste0('analysis/data/derived_data/PITcleanr/TUM_Chinook_', yr, '.rda'))

proc_ch = read_excel(paste0('analysis/data/raw_data/WDFW/PIT_Cleaned_', yr, ' SPCH Data.xlsx')) %>%
  mutate_at(vars(AutoProcStatus:ValidPath),
            list(as.logical)) %>%
  mutate_at(vars(BranchNum, NodeOrder),
            list(as.integer)) %>%
  mutate_at(vars(TrapDate),
            list(lubridate::floor_date),
            unit = "days") %>%
  select(-Group) %>%
  left_join(proc_list$NodeOrder %>%
              select(Node, Group)) %>%
  select(one_of(names(proc_list$ProcCapHist)))

# fix a few double tags
proc_ch %<>%
  left_join(dbl_tag %>%
               select(AltTag = TagID,
                      TagID = TagOther)) %>%
  mutate(TagID = if_else(!is.na(AltTag),
                         AltTag,
                         TagID)) %>%
  select(-AltTag)

proc_ch %>%
  select(TagID, SiteID, Node) %>%
  distinct() %>%
  anti_join(proc_list$ProcCapHist %>%
              select(TagID, Node) %>%
              distinct())

proc_list$ProcCapHist %>%
  select(TagID, SiteID, Node) %>%
  distinct() %>%
  anti_join(proc_ch %>%
              select(TagID, Node) %>%
              distinct())

proc_ch %>%
  filter(Node == 'UWE') %>%
  select(TagID) %>%
  distinct() %>%
  left_join(proc_ch) %>%
  filter(!UserProcStatus) %>%
  select(TagID) %>%
  distinct() %>%
  left_join(proc_ch) %>%
  select(TagID, ObsDate, lastObsDate, SiteID, Node, Direction, matches('Status'))


proc_list$ProcCapHist = proc_ch %>%
  filter(Node != 'UWE')

# re-save some stuff
save(yr, start_date, parent_child, proc_list,
     file = paste0('analysis/data/derived_data/PITcleanr/TUM_Chinook_', yr, '.rda'))


#-----------------------------------------------------------------
# tag summaries
#-----------------------------------------------------------------
bio_df = read_rds('analysis/data/derived_data/Bio_Data_2008_2019.rds')

# Fix UserProcStatus, and summarise tag data
tag_summ = proc_list$ProcCapHist %>%
  filter(UserProcStatus) %>%
  summariseTagData(trap_data = bio_df %>%
                     filter(Year == yr,
                            TagID %in% unique(proc_list$ProcCapHist$TagID)))

# any duplicated tags?
tag_summ %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  as.data.frame()

# where are tags assigned?
# janitor::tabyl(tag_summ, AssignSpawnSite) %>%
janitor::tabyl(tag_summ, AssignSpawnNode) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# which branch are tags assigned to?
tag_summ %>%
  mutate(Branch = fct_explicit_na(Group,
                                 'TUM_bb')) %>%
  janitor::tabyl(Branch) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# preliminary estimate of node efficiency
node_eff = proc_list$ProcCapHist %>%
  filter(AutoProcStatus) %>%
  mutate(UserProcStatus = AutoProcStatus) %>%
  estNodeEff(node_order = proc_list$NodeOrder)

node_eff %>%
  filter(tagsAtNode > 0,
         detEff < 1)

node_eff %>%
  xtabs(~ (!is.na(detEff)) + (detEff_SE > 0), .)

node_eff %>%
  filter(!is.na(detEff),
         detEff_SE > 0)

node_eff %>%
  filter(grepl('^TOP', Node))

#-----------------------------------------------------------------
# examine some of the output
#-----------------------------------------------------------------
proc_ch = proc_list$ProcCapHist

# which tags have "strange" capture histories?
weird_tags = proc_ch %>%
  filter(UserProcStatus == '') %>%
  pull(TagID) %>%
  unique()

length(weird_tags)
length(weird_tags) / n_distinct(proc_ch$TagID)

proc_ch %>%
  filter(TagID %in% weird_tags) %>%
  select(TagID:AutoProcStatus) %>%
  as.data.frame()

proc_ch %>%
  filter(TagID == weird_tags[[1]])


tag_summ %>%
  mutate(Population = Group,
         Population = if_else(AssignSpawnSite == 'LWC' | grepl('ROZ', TagPath),
                              'Upper Yakima',
                              Population),
         Population = if_else(AssignSpawnSite %in% c('LNR', 'AH1'),
                              'Naches',
                              Population),
         Population = if_else(is.na(Group),
                              'PRO_bb',
                              Population)) %>%
  filter(Population %in% c('Status', 'Naches', 'Toppenish', 'Upper Yakima')) %>%
  ggplot(aes(x = PassDate,
             color = Population,
             fill = Population)) +
  # geom_density(alpha = 0.2) +
  geom_histogram() +
  facet_wrap(~ Population) +
  theme_bw() +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1')

# ggsave('outgoing/figures/RunTiming_2019.pdf',
#        width = 6,
#        height = 6)

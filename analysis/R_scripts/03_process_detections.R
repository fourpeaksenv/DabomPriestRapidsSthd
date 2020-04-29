# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 4/27/20
# Last Modified: 4/27/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(lubridate)
library(readxl)
library(magrittr)

#-----------------------------------------------------------------
# load configuration and site_df data
load('analysis/data/derived_data/site_config.rda')

# which spawn year are we dealing with?
yr = 2016

# start date is June 1 of the previous year
start_date = paste0(yr - 1, '0601')

# build parent-child table
parent_child = createParentChildDf(site_df,
                                   configuration,
                                   startDate = start_date)

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
observations = read_csv(paste0('analysis/data/raw_data/PTAGIS/UC_Sthd_', yr, '_CTH.csv'))

# process those observations with PITcleanr, using Tumwater-specific function
proc_list = processCapHist_PRD(startDate = start_date,
                               configuration = configuration,
                               parent_child = parent_child,
                               observations = observations,
                               last_obs_date = format(ymd(start_date) + years(1) + months(1), "%Y%m%d"),
                               truncate = T,
                               site_df = site_df,
                               step_num = 1,
                               save_file = T,
                               file_name = paste0('outgoing/PITcleanr/UC_Steelhead_', yr, '.xlsx'))

# create node order and list of nodes within sevaral population groups
node_order = createNodeOrder(proc_list$ValidPaths,
                             configuration) %>%
  left_join(stack(site_list) %>%
              tbl_df() %>%
              select(Group = ind,
                     NodeSite = values) %>%
              mutate(BranchNum = as.integer(Group))) %>%
  distinct()

proc_list$NodeOrder = node_order

# save some stuff
save(yr, start_date, parent_child, proc_list,
     file = paste0('analysis/data/derived_data/PITcleanr/UC_Steelhead_', yr, '.rda'))

#-------------------------------------------
# NEXT STEPS
#-------------------------------------------
# open that Excel file, and filter on the column UserProcStatus, looking for blanks. Fill in each row with TRUE or FALSE, depending on whether that observation should be kept or not. The column AutoProcStatus provides a suggestion, but the biologist's best expert judgement should be used.

#-------------------------------------------
# After receiving cleaned up file back...
#-------------------------------------------
load(paste0('analysis/data/derived_data/PITcleanr/UC_Steelhead_', yr, '.rda'))

proc_ch = read_excel(paste0('analysis/data/raw_data/WDFW/UC_Steelhead_', yr, '.xlsx')) %>%
  mutate_at(vars(AutoProcStatus:ValidPath),
            list(as.logical)) %>%
  mutate_at(vars(BranchNum, NodeOrder),
            list(as.integer)) %>%
  mutate_at(vars(TrapDate),
            list(lubridate::ymd)) %>%
  mutate_at(vars(ObsDate, lastObsDate),
            list(lubridate::ymd_hms)) %>%
  mutate_at(vars(TrapDate),
            list(lubridate::floor_date),
            unit = "days") %>%
  select(-Group) %>%
  left_join(proc_list$NodeOrder %>%
              select(Node, Group)) %>%
  select(one_of(names(proc_list$ProcCapHist)))


proc_ch %>%
  select(TagID, SiteID, Node) %>%
  distinct() %>%
  anti_join(proc_list$ProcCapHist %>%
              select(TagID, Node) %>%
              distinct()) %>%
  xtabs(~ SiteID, .)

proc_list$ProcCapHist %>%
  select(TagID, SiteID, Node) %>%
  distinct() %>%
  anti_join(proc_ch %>%
              select(TagID, Node) %>%
              distinct()) %>%
  xtabs(~ SiteID, .)

# overwrite previously saved capture histories
# proc_list$ProcCapHist = proc_ch

# re-save some stuff
save(yr, start_date, parent_child, proc_list,
     file = paste0('analysis/data/derived_data/PITcleanr/UC_Steelhead_', yr, '.rda'))


#-----------------------------------------------------------------
# tag summaries
#-----------------------------------------------------------------
file_nms = list.files('analysis/data/derived_data')
bio_nm = file_nms[grepl('Bio', file_nms) & grepl('.rds$', file_nms)]

bio_df = read_rds(paste0('analysis/data/derived_data/', bio_nm))

# Fix UserProcStatus, and summarise tag data
tag_summ = proc_list$ProcCapHist %>%
  filter(UserProcStatus) %>%
  summariseTagData(trap_data = bio_df %>%
                     filter(Year == yr,
                            TagID %in% unique(proc_list$ProcCapHist$TagID)))

# any duplicated tags?
tag_summ %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  arrange(TagID, TrapDate) %>%
  as.data.frame()

tag_summ %<>%
  group_by(TagID) %>%
  filter(TrapDate == min(TrapDate))

# where are tags assigned?
janitor::tabyl(tag_summ, AssignSpawnSite) %>%
# janitor::tabyl(tag_summ, AssignSpawnNode) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# which branch are tags assigned to?
tag_summ %>%
  mutate(Branch = fct_explicit_na(Group,
                                 'PRA_bb')) %>%
  janitor::tabyl(Branch, Origin) %>%
  # janitor::adorn_totals(where = c("row")) %>%
  # janitor::adorn_percentages(denominator = "col") %>%
  # janitor::adorn_pct_formatting()
  janitor::adorn_totals(where = c("row", "col"))

# preliminary estimate of node efficiency
node_eff = proc_list$ProcCapHist %>%
  filter(UserProcStatus) %>%
  estNodeEff(node_order = proc_list$NodeOrder)

node_eff %>%
  filter(tagsAtNode > 0,
         detEff < 1)

node_eff %>%
  filter(!is.na(detEff),
         detEff_SE > 0)

node_eff %>%
  filter(grepl('^LWE', Node))

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

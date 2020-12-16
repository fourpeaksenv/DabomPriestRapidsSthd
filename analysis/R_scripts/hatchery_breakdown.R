# Author: Kevin See
# Purpose: break down hatchery estimates into different groups
# Created: 12/8/20
# Last Modified: 12/8/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(coda)
library(PITcleanr)
library(DABOM)
library(readxl)

#-----------------------------------------------------------------
# set species
spp = "Steelhead"
# set year
yr = 2019

# load JAGS MCMC results
load(paste0("analysis/data/derived_data/model_fits/PRA_", spp, "_", yr,'_DABOM.rda'))

# estimate final spawning location
tag_summ = summariseTagData(proc_list$proc_ch %>%
                              select(-Group) %>%
                              left_join(proc_list$NodeOrder %>%
                                          select(Node, Group)),
                            trap_data = bio_df %>%
                              filter(TagID %in% unique(proc_list$proc_ch$TagID)) %>%
                              group_by(TagID) %>%
                              slice(1)) %>%
  mutate(Group = fct_explicit_na(Group))


# any duplicate tags?
tag_summ %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  arrange(TagID, TrapDate)

# get valid paths by site (not node)
valid_paths = parent_child %>%
  mutate(across(c(ParentNode,
                  ChildNode),
                str_remove,
                pattern = "A0"),
         across(c(ParentNode,
                  ChildNode),
                str_remove,
                pattern = "B0")) %>%
  filter(ParentNode != ChildNode) %>%
  getValidPaths(root_site = "PRA") %>%
  rename(Site = Node)



tag_summ %>%
  mutate(CWT = if_else(is.na(CWT), F,
                        if_else(CWT == "SN", T, NA)),
         AdClip = if_else(is.na(AdClip), F,
                          if_else(AdClip == "AD", T, NA))) %>%
  mutate(Site = str_remove(AssignSpawnNode, "B0$"),
         Site = str_remove(Site, "A0$")) %>%
  left_join(valid_paths %>%
              mutate(path_df = str_split(Path, " "),
                     path_df = map(path_df,
                                   .f = function(x) {
                                     tibble(path_sites = x)
                                   })) %>%
              select(-Path) %>%
              unnest(path_df)) %>%
  group_by(path_sites, AdClip, CWT) %>%
  summarise(n_tags = n_distinct(TagID),
            .groups = "drop") %>%
  group_by(path_sites) %>%
  mutate(prop = n_tags / sum(n_tags))



# get DABOM estimates of hatchery escapement
hatch_est = read_excel(paste0('outgoing/estimates/PRA_', spp, '_', yr, '_20200604.xlsx'),
                       sheet = 'All Escapement') %>%
  filter(Origin == "Hatchery",
         grepl("^past_", param),
         estimate > 0) %>%
  mutate(node = str_remove(param, "^past_"))

# Author: Kevin See
# Purpose: break down hatchery estimates into different groups
# Created: 12/8/20
# Last Modified: 12/17/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(coda)
library(PITcleanr)
library(DABOM)
library(STADEM)
library(readxl)
library(msm)

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
sum(duplicated(tag_summ$TagID))
tag_summ %>%
  filter(TagID %in% TagID[duplicated(TagID)]) %>%
  arrange(TagID, TrapDate)

# get valid paths by site (not node)
valid_paths = parent_child %>%
  mutate(across(c(ParentNode,
                  ChildNode),
                str_remove,
                pattern = "A0$"),
         across(c(ParentNode,
                  ChildNode),
                str_remove,
                pattern = "B0$")) %>%
  # change one node back
  mutate(across(c(ParentNode,
                  ChildNode),
                  ~ recode(., "S" = "SA0"))) %>%
  filter(ParentNode != ChildNode) %>%
  getValidPaths(root_site = "PRA") %>%
  rename(Site = Node)

# proportions of each type of fish (Ad-clip, CWT combinations) past each site
tag_prop = tag_summ %>%
  filter(Origin == "H") %>%
  mutate(CWT = if_else(is.na(CWT), F,
                        if_else(CWT %in% c("SN", "BD"),
                                T, NA)),
         AdClip = if_else(is.na(AdClip), F,
                          if_else(AdClip == "AD", T, NA))) %>%
  mutate(Site = str_remove(AssignSpawnNode, "B0$"),
         Site = str_remove(Site, "A0$"),
         Site = recode(Site,
                       "S" = "SA0")) %>%
  left_join(valid_paths %>%
              mutate(path_df = str_split(Path, " "),
                     path_df = map(path_df,
                                   .f = function(x) {
                                     tibble(path_sites = x)
                                   })) %>%
              select(-Path) %>%
              unnest(path_df)) %>%
  mutate(path_sites = if_else(is.na(path_sites) & AssignSpawnNode == "PRA",
                              "PRA",
                              path_sites)) %>%
  group_by(path_sites, AdClip, CWT) %>%
  summarise(n_tags = n_distinct(TagID),
            .groups = "drop") %>%
  right_join(., expand(., path_sites, AdClip, CWT)) %>%
  arrange(path_sites, AdClip, CWT) %>%
  mutate(across(n_tags,
                replace_na,
                replace = 0)) %>%
  group_by(path_sites) %>%
  mutate(tot_tags = sum(n_tags),
         prop = n_tags / tot_tags,
         # using normal approximation
         prop_se = sqrt((prop * (1 - prop))/tot_tags)) %>%
  ungroup() %>%
  rename(Site = path_sites)

# need to figure out multinomial proportion standard error

# get DABOM estimates of hatchery escapement
# movement probabilities
trans_df = compileTransProbs_PRA(dabom_mod) %>%
  filter(Origin == "Hatchery")

# get escapement past Priest by origin
start_date = paste0(yr-1, '0601')
end_date = paste0(yr, '0531')

org_escape = queryPITtagData(damPIT = 'PRA',
                             spp = spp,
                             start_date = start_date,
                             end_date = end_date) %>%
  mutate(SpawnYear = yr,
         TagIDAscentCount = ifelse(is.na(TagIDAscentCount),
                                   0, TagIDAscentCount),
         ReAscent = ifelse(TagIDAscentCount > 1, T, F)) %>%
  group_by(Species, SpawnYear, Date) %>%
  summarise(tot_tags = n_distinct(TagID),
            reascent_tags = n_distinct(TagID[ReAscent]),
            .groups = "drop") %>%
  group_by(Species, SpawnYear) %>%
  summarise(across(matches('tags'),
                   sum,
                   na.rm = T),
            .groups = "drop") %>%
  mutate(reascRate = reascent_tags / tot_tags,
         reascRateSE = sqrt(reascRate * (1 - reascRate) / tot_tags)) %>%
  bind_cols(getWindowCounts(dam = 'PRD',
                            spp = spp,
                            start_date = paste0(yr-1, '0601'),
                            end_date = paste0(yr, '0531')) %>%
              summarise_at(vars(win_cnt),
                           list(sum)) %>%
              select(tot_win_cnt = win_cnt)) %>%
  mutate(adjWinCnt = tot_win_cnt * (1 - reascRate),
         adjWinCntSE = tot_win_cnt * reascRateSE) %>%
  bind_cols(bio_df %>%
              group_by(Origin) %>%
              summarise(nTags = n_distinct(TagID)) %>%
              pivot_wider(names_from = "Origin",
                          values_from = "nTags",
                          values_fill = list(nTags = as.integer(0))) %>%
              ungroup() %>%
              mutate(propW = W / (W + H),
                     propH = 1 - propW,
                     propOrgSE = sqrt((propW * (1 - propW)) / (W + H)))) %>%
  mutate(Hescp = propH * adjWinCnt,
         HescpSE = msm::deltamethod(~ x1 * x2,
                               mean = c(propH, adjWinCnt),
                               cov = diag(c(propOrgSE, adjWinCntSE)^2))) %>%
  mutate(Wescp = propW * adjWinCnt,
         WescpSE = msm::deltamethod(~ x1 * x2,
                               mean = c(propW, adjWinCnt),
                               cov = diag(c(propOrgSE, adjWinCntSE)^2))) %>%
  select(Species, SpawnYear, matches('escp')) %>%
  pivot_longer(-(Species:SpawnYear),
               names_to = "var",
               values_to = "value") %>%
  mutate(Origin = if_else(grepl('^H', var),
                          'Hatchery',
                          'Natural'),
         param = if_else(grepl('SE$', var),
                         'tot_escp_se',
                         'tot_escp')) %>%
  select(-var) %>%
  pivot_wider(names_from = "param",
              values_from = "value")

set.seed(3)
escp_post = org_escape %>%
  filter(Origin == "Hatchery") %>%
  crossing(chain = unique(trans_df$chain)) %>%
  mutate(samps = map2(tot_escp,
                      tot_escp_se,
                      .f = function(x, y) {
                        tibble(escp = rnorm(max(trans_df$iter),
                                            x, y)) %>%
                          mutate(iter = 1:n())
                      })) %>%
  select(Species, Origin, chain, samps) %>%
  unnest(cols = samps) %>%
  inner_join(trans_df %>%
               filter(grepl('^past', param))) %>%
  mutate(value = value * escp) %>%
  mutate(Site = str_remove(param, "past_"))

set.seed(6)
prop_samps = tag_prop %>%
  crossing(chain = unique(trans_df$chain)) %>%
  mutate(across(starts_with("prop"),
                ~ if_else(. == 0, 1e-10, .))) %>%
  mutate(prop_alpha = ((1 - prop) / prop_se^2 - 1 / prop) * prop^2,
         prop_alpha = if_else(prop_alpha < 0, 0.01, prop_alpha),
         prop_beta = prop_alpha * (1 / prop - 1)) %>%
  mutate(samps = map2(prop_alpha,
                      prop_beta,
                      .f = function(x, y) {
                        tibble(prop = rbeta(max(trans_df$iter),
                                            x, y)) %>%
                          mutate(iter = 1:n())
                      })) %>%
  select(Site:CWT,
         n_tags, tot_tags,
         # prop_obs = prop,
         # prop_se,
         chain, samps) %>%
  # filter(prop_obs > 0,
  #        prop_obs < 1) %>%
  unnest(cols = samps)


all_post = escp_post %>%
  inner_join(prop_samps) %>%
  mutate(n_fish = value * prop) %>%
  mutate(across(n_fish,
                round,
                digits = 2))

mark_grps = all_post %>%
  group_by(Species,
           Origin,
           Site,
           AdClip,
           CWT,
           n_tags, tot_tags) %>%
  summarise(across(n_fish,
                   tibble::lst(mean,
                               median,
                               mode = estMode,
                               se = sd,
                               skew = moments::skewness,
                               kurtosis = moments::kurtosis)),
            .groups = "drop") %>%
  full_join(all_post %>%
              group_by(Species,
                       Origin,
                       Site,
                       AdClip,
                       CWT) %>%
              summarise(lowerCI = coda::HPDinterval(coda::as.mcmc(n_fish))[,1],
                        upperCI = coda::HPDinterval(coda::as.mcmc(n_fish))[,2],
                        .groups = "drop"))

names(mark_grps) = str_remove(names(mark_grps), "^n_fish_")

mark_grps %>%
  filter(n_tags > 0)

mark_grps %>%
  full_join(hatch_est %>%
              select(Origin, Site,
                     tot_escp = estimate,
                     tot_se = se)) %>%
  tail()





hatch_est = read_excel(paste0('outgoing/estimates/PRA_', spp, '_', yr, '_20200604.xlsx'),
                       sheet = 'All Escapement') %>%
  filter(Origin == "Hatchery",
         !grepl('_bb$', param),
         # grepl("^past_", param),
         estimate > 0) %>%
  mutate(Site = str_remove(param, "^past_"))


hatch_type = hatch_est %>%
  select(Origin, Site,
         tot_escp = estimate,
         tot_se = se) %>%
  full_join(tag_prop %>%
              filter(Site != "PRA") %>%
              select(Site:CWT,
                     tot_tags,
                     n_tags,
                     prop)) %>%
  mutate(N_hat = tot_escp * prop) %>%
  full_join(mark_grps) %>%
  # filter(n_tags > 0) %>%
  mutate(across(c(N_hat,
                  mean:mode),
                round))

hatch_type %>%
  mutate(across(c(AdClip, CWT),
                as.character)) %>%
  mutate(AdClip = recode(AdClip,
                         "TRUE" = "AD",
                         "FALSE" = "AI"),
         CWT = recode(CWT,
                      "TRUE" = "CWT",
                      "FALSE" = "noCWT")) %>%
  pivot_wider(id_cols = c(Species, Origin:tot_se),
              names_from = c("AdClip", "CWT"),
              values_from = "mean",
              values_fill = 0) %>%
  filter(is.na(Species))

# combine total hatchery escapement with proportions of CWT / Ad Clip combos
hatch_type = hatch_est %>%
  left_join(tag_prop) %>%
  rowwise() %>%
  mutate(escp = estimate * prop,
         escp_se = msm::deltamethod(~ x1 * x2,
                                    mean = c(estimate, prop),
                                    cov = diag(c(se, prop_se)^2))) %>%
  ungroup() %>%
  select(Origin, Site,
         AdClip,
         CWT,
         tot_tags,
         n_tags,
         tot_escp = estimate,
         tot_se = se,
         prop:escp_se)

hatch_type %>%
  pivot_wider(id_cols = c(Origin:tot_se),
              names_from = c("AdClip", "CWT"),
              values_from = "escp") %>%
  filter(is.na(tot_tags) | tot_tags == 0)

hatch_est$Site[!hatch_est$Site %in% hatch_type$Site]

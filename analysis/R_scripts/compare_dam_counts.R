# Author: Kevin See
# Purpose: recalculate estimates at Priest based on other dam counts
# Created: 1/25/2022
# Last Modified: 1/25/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(STADEM)
library(here)
library(DABOM)
library(lubridate)
library(janitor)
library(msm)
library(magrittr)

theme_set(theme_bw())

#-----------------------------------------------------------------
# what years to examine?
all_yrs = 2011:2021

# set species
spp = "Steelhead"

#-----------------------------------------------------------------
# get the total dam counts from various dams
dam_cnts = tibble(dam = c("PriestRapids",
                          "RockIsland",
                          "RockyReach",
                          "Wells",
                          "Tumwater"),
                  dam_code = c("PRD",
                               "RIS",
                               "RRH",
                               "WEL",
                               "TUM")) %>%
  crossing(year = all_yrs) %>%
  rowwise() %>%
  mutate(win_cnt = map2_dbl(dam_code,
                            year,
                           .f = function(x, y) {
                             STADEM::getWindowCounts(dam = x,
                                                     spp = spp,
                                                     start_date = paste0(y-1, '0601'),
                                                     end_date = paste0(y, '0531')) %>%
                               summarise_at(vars(win_cnt),
                                            list(sum),
                                            na.rm = T) %>%
                               pull(win_cnt)
                           })) %>%
  ungroup() %>%
  select(year, everything())

# got these directly from Ben for each spawn year (June 15 - June 14, which is how it's definted at Tumwater)
tum_cnts = tibble(year = 2013:2021,
                  dam = "Tumwater",
                  dam_code = "TUM",
                  win_cnt = c(2446,
                              1186,
                              1751,
                              1405,
                              554,
                              621,
                              390,
                              578,
                              776))

dam_cnts %<>%
  left_join(tum_cnts %>%
              rename(tum_cnt = win_cnt)) %>%
  mutate(win_cnt = if_else(!is.na(tum_cnt),
                           tum_cnt,
                           win_cnt)) %>%
  select(-tum_cnt)

#-----------------------------------------------------------------
# loop over all years
for(yr in all_yrs) {

  # set up a tibble to capture results
  if(yr == min(all_yrs)) {
    pit_move_all = NULL
  }

  if(spp == "Steelhead") {
    start_date = paste0(yr-1, '0601')
    end_date = paste0(yr, '0531')
  }

  # load compressed detections and biological data
  load(here('analysis/data/derived_data/PITcleanr',
            paste0('UC_', spp, '_', yr, '.rda')))


  # load JAGS MCMC results
  load(here("analysis/data/derived_data/model_fits",
            paste0('PRA_DABOM_', spp, '_', yr,'.rda')))

  # summarize detection probabilities
  detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                     filter_ch = filter_obs)
  # add origin to prepped capture histories
  pit_obs = prepped_ch %>%
    left_join(bio_df %>%
                select(tag_code,
                       origin)) %>%
    select(tag_code, origin,
           everything())

  # compile and calculate Priest equivalents
  pit_move_df = tibble(dam = c("PriestRapids",
                               "RockIsland",
                               "RockyReach",
                               "Wells",
                               "Tumwater"),
                       dam_code = c("PRD",
                                    "RIS",
                                    "RRH",
                                    "WEL",
                                    "TUM"),
                       pit_code = c("PRA",
                                    "RIA",
                                    "RRF",
                                    "WEA",
                                    "TUM")) %>%
    add_column(year = yr,
               .before = 0) %>%
    crossing(origin = c("W", "H")) %>%
    left_join(pit_obs %>%
                group_by(origin) %>%
                summarize(n_tags = n_distinct(tag_code))) %>%
    mutate(n_obs = map2_int(pit_code,
                            origin,
                            .f = function(x, y) {
                              pit_obs %>%
                                filter(origin == y) %>%
                                summarize(n_obs = n_distinct(tag_code[node == x])) %>%
                                pull(n_obs)
                            })) %>%
    mutate(n_upstrm = map2_int(pit_code,
                               origin,
                               .f = function(x, y) {
                                 pit_obs %>%
                                   filter(origin == y) %>%
                                   summarize(n_path = n_distinct(tag_code[str_detect(path, paste0(" ", x))])) %>%
                                   pull(n_path)
                               }),
           n_upstrm = if_else(dam == "PriestRapids",
                              n_obs,
                              n_upstrm)) %>%
    # generate proportion hatchery/wild based on upstream tags
    group_by(dam) %>%
    mutate(prop_org = n_upstrm / sum(n_upstrm),
           prop_org_se = sqrt((prop_org * (1 - prop_org)) / sum(n_upstrm))) %>%
    # generate proportion hatchery/wild based on tags detected at dam
    # mutate(prop_org = n_obs / sum(n_obs),
    #        prop_org_se = sqrt((prop_org * (1 - prop_org)) / sum(n_obs))) %>%
    ungroup() %>%
    # mutate(prop_obs = n_obs / n_tags,
    #        prop_upstrm = n_upstrm / n_tags) %>%
    left_join(detect_summ %>%
                select(node, mean, sd),
              by = c("pit_code" = "node")) %>%
    mutate(est_tags = n_obs / mean,
           est_tags_se = n_obs * sd / mean^2,
           trans_est = est_tags / n_tags,
           trans_se = est_tags_se / n_tags) %>%
    mutate(trans_est = if_else(dam == "PriestRapids",
                               1,
                               trans_est),
           trans_se = if_else(dam == "PriestRapids",
                              0,
                              trans_se)) %>%
    left_join(dam_cnts) %>%
    rowwise() %>%
    mutate(priest_cnt = win_cnt * prop_org / trans_est,
           priest_cnt_se = msm::deltamethod(~ x1 * x2 / x3,
                                            mean = c(win_cnt,
                                                     prop_org,
                                                     trans_est),
                                            cov = diag(c(0,
                                                         prop_org_se,
                                                         trans_se)^2))) %>%
    ungroup() %>%
    # add re-ascenstion rate at Priest
    left_join(queryPITtagData(damPIT = 'PRA',
                              spp = "Steelhead",
                              start_date = start_date,
                              end_date = end_date) %>%
                filter(!str_detect(TagId, "000.0")) %>%
                mutate(SpawnYear = yr) %>%
                mutate(across(TagIdAscentCount,
                              tidyr::replace_na,
                              0)) %>%
                mutate(ReAscent = ifelse(TagIdAscentCount > 1, T, F)) %>%
                mutate(origin = fct_recode(RearType,
                                           "W" = "U")) %>%
                group_by(Species, SpawnYear, Date, origin) %>%
                summarise(tot_tags = n_distinct(TagId),
                          reascent_tags = n_distinct(TagId[ReAscent]),
                          .groups = "drop") %>%
                group_by(Species, SpawnYear, origin) %>%
                summarise(across(matches('tags'),
                                 sum,
                                 na.rm = T),
                          .groups = "drop") %>%
                mutate(reasc_rate = reascent_tags / tot_tags,
                       reasc_rate_se = sqrt(reasc_rate * (1 - reasc_rate) / tot_tags)) %>%
                select(origin, starts_with("reasc_rate"))) %>%
    rowwise() %>%
    mutate(tot_escp = priest_cnt * (1 - reasc_rate),
           tot_escp_se = msm::deltamethod(~ x1 * (1 - x2),
                                          mean = c(priest_cnt,
                                                   reasc_rate),
                                          cov = diag(c(priest_cnt_se,
                                                       reasc_rate_se)^2))) %>%
    ungroup() %>%
    mutate(dam = fct_relevel(dam,
                             "Tumwater",
                             after = Inf))

  if(is.null(pit_move_all)) {
    pit_move_all = pit_move_df
  } else {
    pit_move_all <- pit_move_all %>%
      bind_rows(pit_move_df)
  }

  rm(detect_summ,
     pit_obs,
     pit_move_df,
     start_date, end_date)
}

pit_move_all %<>%
    mutate(priest_cnt_se = if_else(dam == "PriestRapids",
                                   NA_real_,
                                   priest_cnt_se))


#-----------------------------------------------------------------
# make a density plot
set.seed(5)
dens_comp_p = pit_move_all %>%
  select(year, dam, origin, win_cnt,
         matches('priest_cnt')) %>%
  group_by(year, dam, win_cnt) %>%
  summarise(priest_cnt = sum(priest_cnt),
            priest_cnt_se = sqrt(sum(priest_cnt_se^2)),
            .groups = "drop") %>%
  mutate(priest_cnt_se = if_else(dam == "PriestRapids",
                                 NA_real_,
                                 priest_cnt_se)) %>%
  # mutate(across(win_cnt:priest_cnt_se,
  #               ~ . / 1000)) %>%
  filter(dam != "PriestRapids") %>%
  mutate(samps = map2(priest_cnt,
                      priest_cnt_se,
                      .f = function(x, y) rnorm(10000, x, y))) %>%
  unnest(samps) %>%
  ggplot(aes(x = samps,
             color = dam,
             fill = dam)) +
  scale_color_brewer(palette = "Set2",
                     name = "Dam Count\nSource") +
  scale_fill_brewer(palette = "Set2",
                    name = "Dam Count\nSource") +
  geom_density(alpha = 0.4) +
  geom_vline(data = dam_cnts %>%
               filter(dam == "PriestRapids"),
             aes(xintercept = win_cnt),
             linetype = 2,
             lwd = 1) +
  facet_wrap(~ year,
             scales = "free") +
  scale_x_continuous(labels = function(x) x / 1000) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = "Equivalent of Total Counts\nat Priest Rapids Dam (x1,000)")
dens_comp_p

xy_p = pit_move_all %>%
  select(year, dam, origin, win_cnt,
         matches('priest_cnt')) %>%
  group_by(year, dam, win_cnt) %>%
  summarise(priest_cnt = sum(priest_cnt),
            priest_cnt_se = sqrt(sum(priest_cnt_se^2)),
            .groups = "drop") %>%
  mutate(priest_cnt_se = if_else(dam == "PriestRapids",
                                 NA_real_,
                                 priest_cnt_se)) %>%
  ggplot(aes(x = dam,
             y = priest_cnt,
             color = dam)) +
  geom_errorbar(aes(ymax = qnorm(0.975, priest_cnt, priest_cnt_se),
                    ymin = qnorm(0.025, priest_cnt, priest_cnt_se)),
                width = 0.3) +
  geom_point(size = 4) +
  facet_wrap(~ year,
             scales = "free_y") +
  scale_color_brewer(palette = "Set1",
                     name = "Dam Count Source") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  scale_y_continuous(labels = function(x) x / 1000) +
  labs(x = "Dam Count Source",
       y = "Estimated Counts at Priest (x1,000)")
xy_p

rel_xy_p = pit_move_all %>%
  select(year, dam, origin, win_cnt,
         matches('priest_cnt')) %>%
  group_by(year, dam, win_cnt) %>%
  summarise(priest_cnt = sum(priest_cnt),
            priest_cnt_se = sqrt(sum(priest_cnt_se^2)),
            .groups = "drop") %>%
  mutate(priest_cnt_se = if_else(dam == "PriestRapids",
                                 NA_real_,
                                 priest_cnt_se)) %>%
  group_by(year) %>%
  mutate(across(starts_with("priest_cnt"),
                ~ . / win_cnt[dam == "PriestRapids"]),
         across(priest_cnt,
                ~ . - 1)) %>%
  ggplot(aes(x = dam,
             y = priest_cnt,
             color = dam)) +
  geom_errorbar(aes(ymax = qnorm(0.975, priest_cnt, priest_cnt_se),
                    ymin = qnorm(0.025, priest_cnt, priest_cnt_se)),
                width = 0.3) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0,
             linetype = 2) +
  facet_wrap(~ year,
             # scales = "fixed") +
             scales = "free_y") +
  scale_color_brewer(palette = "Set1",
                     name = "Dam Count Source") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  scale_y_continuous(labels = function(x) x * 100) +
  labs(x = "Dam Count Source",
       y = "Relative Difference\nCompared to Priest (%)")
rel_xy_p

ggsave(here("outgoing/other/dam_cnt_comp_density.pdf"),
       dens_comp_p,
       width = 8,
       height = 6)

ggsave(here("outgoing/other/dam_cnt_comp_xy.pdf"),
       xy_p,
       width = 8,
       height = 6)

ggsave(here("outgoing/other/dam_cnt_comp_rel.pdf"),
       rel_xy_p,
       width = 8,
       height = 6)

d_wd = 0.5
# pit_move_all %>%
#   select(year, dam, origin, win_cnt,
#          matches('priest_cnt')) %>%
#   group_by(year, dam, win_cnt) %>%
#   summarise(priest_cnt = sum(priest_cnt),
#             priest_cnt_se = sqrt(sum(priest_cnt_se^2)),
#             .groups = "drop") %>%
#   mutate(origin = "All") %>%
#   mutate(priest_cnt_se = if_else(dam == "PriestRapids",
#                                  NA_real_,
#                                  priest_cnt_se)) %>%
#   bind_rows(pit_move_all) %>%
#   arrange(year, origin, dam) %>%
pit_move_all %>%
  ggplot(aes(x = as_factor(year),
             y = priest_cnt,
             color = dam)) +
  geom_errorbar(aes(ymax = qnorm(0.975, priest_cnt, priest_cnt_se),
                    ymin = qnorm(0.025, priest_cnt, priest_cnt_se)),
                width = 0.1,
                position = position_dodge(width = d_wd)) +
  geom_point(size = 2,
             position = position_dodge(width = d_wd)) +
  facet_wrap(~ origin,
             # scales = "fixed") +
             scales = "free_y") +
  scale_color_brewer(palette = "Set1",
                     name = "Dam Count Source") +
  theme(axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  labs(x = "Dam Count Source",
       y = "Estimated Counts at Priest")

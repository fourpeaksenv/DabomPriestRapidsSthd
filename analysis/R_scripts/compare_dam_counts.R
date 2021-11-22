# Author: Kevin See
# Purpose: recalculate estimates at Priest based on other dam counts
# Created: 9/22/21
# Last Modified: 9/22/21
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(STADEM)
library(here)
library(DABOM)

theme_set(theme_bw())

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data',
          'site_config.rda'))

#-----------------------------------------------------------------
# set a year
yr = 2018

start_date = paste0(yr-1, '0601')
end_date = paste0(yr, '1231')

# query dam counts at Priest, Rock Island, Rocky Reach and Wells
dam_cnts = list(PriestRapids = "PRD",
                RockIsland = "RIS",
                RockyReach = "RRH",
                Wells = "WEL",
                Tumwater = "TUM") %>%
  map_df(.id = "dam",
         .f = function(x) {
           getWindowCounts(dam = x,
                           spp = "Steelhead",
                           start_date = start_date,
                           end_date = end_date) %>%
             summarise_at(vars(win_cnt),
                          list(sum),
                          na.rm = T)
         })

dam_cnts


# load JAGS MCMC results
load(here("analysis/data/derived_data/model_fits",
          paste0('PRA_DABOM_Steelhead_', yr,'.rda')))

# cdam_cntsompile all movement probabilities, and multiply them appropriately
trans_df = compileTransProbs_PRA(dabom_mod,
                                 parent_child) %>%
  filter(param %in% c("RIA",
                      "RRF",
                      "WEA",
                      "TUM"),
         origin == 2) %>%
  mutate(dam = recode(param,
                      "RIA" = "RockIsland",
                      "RRF" = "RockyReach",
                      "WEA" = "Wells",
                      "TUM" = "Tumwater"))

priest_post = trans_df %>%
  left_join(dam_cnts) %>%
  mutate(priest_cnt = win_cnt / value)

priest_summ = priest_post %>%
  group_by(dam) %>%
  summarise(mean = mean(priest_cnt),
            median = median(priest_cnt),
            mode = estMode(priest_cnt),
            sd = sd(priest_cnt),
            skew = moments::skewness(priest_cnt),
            kurtosis = moments::kurtosis(priest_cnt),
            lowerCI = coda::HPDinterval(coda::as.mcmc(priest_cnt))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(priest_cnt))[,2],
            .groups = 'drop')

#-----------------------------------------------------------------
# create a plot
priest_p = priest_post %>%
  ggplot(aes(x = priest_cnt,
             color = dam,
             fill = dam)) +
  geom_density(alpha = 0.2) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_vline(xintercept = dam_cnts %>%
               filter(dam == "PriestRapids") %>%
               pull(win_cnt),
             lwd = 2,
             linetype = 2) +
  labs(x = "Estimate",
       color = "Based on Counts at:",
       fill = "Based on Counts at:",
       title = paste("Steelhead Crossing Priest Rapids in", yr))
priest_p

ggsave(here("analysis/figures",
            paste0("PriestEstimates_", yr, ".pdf")),
       width = 6,
       height = 6)

#-----------------------------------------------------------------
# compare dam counts to model estimates
save_list = read_rds(here('analysis/data/derived_data/results',
                          paste0("PRA_Steelhead_estimates_", yr, ".rds")))
save_list$`All Escapement` %>%
  filter(location %in% c("RIA", "RRF", "WEA", "TUM")) %>%
  left_join(dam_cnts %>%
              filter(dam != "PriestRapids") %>%
              mutate(location = recode(dam,
                                       "RockIsland" = "RIA",
                                       "RockyReach" = "RRF",
                                       "Wells" = "WEA",
                                       "Tumwater" = "TUM")))


#-----------------------------------------------------------------
# select a series of years
#-----------------------------------------------------------------
comp_df = tibble(year = 2011:2020) %>%
  mutate(dam_cnts = map(year,
                        .f = function(yr) {

                          start_date = paste0(yr-1, '0601')
                          end_date = paste0(yr, '1231')

                          # query dam counts at Priest, Rock Island, Rocky Reach and Wells
                          dam_cnts = list(PriestRapids = "PRD",
                                          Wanapum = "WAN",
                                          RockIsland = "RIS",
                                          RockyReach = "RRH",
                                          Wells = "WEL",
                                          Tumwater = "TUM") %>%
                            map_df(.id = "dam",
                                   .f = function(x) {
                                     getWindowCounts(dam = x,
                                                     spp = "Steelhead",
                                                     start_date = start_date,
                                                     end_date = end_date) %>%
                                       summarise_at(vars(win_cnt),
                                                    list(sum),
                                                    na.rm = T)
                                   })

                          return(dam_cnts)
                        })) %>%
  mutate(priest_cnts = map_dbl(dam_cnts,
                               .f = function(x) {
                                 x %>%
                                   filter(dam == "PriestRapids") %>%
                                   pull(win_cnt)
                               })) %>%
  mutate(dabom_est = map2(year, dam_cnts,
                         .f = possibly(function(yr, dam_cnts) {

                           # load JAGS MCMC results
                           load(here("analysis/data/derived_data/model_fits",
                                     paste0('PRA_DABOM_Steelhead_', yr,'.rda')))

                           # cdam_cntsompile all movement probabilities, and multiply them appropriately
                           trans_df = compileTransProbs_PRA(dabom_mod,
                                                            parent_child) %>%
                             filter(param %in% c("RIA",
                                                 "RRF",
                                                 "WEA",
                                                 "TUM"),
                                    origin == 2) %>%
                             mutate(dam = recode(param,
                                                 "RIA" = "RockIsland",
                                                 "RRF" = "RockyReach",
                                                 "WEA" = "Wells",
                                                 "TUM" = "Tumwater"))

                           priest_post = trans_df %>%
                             left_join(dam_cnts) %>%
                             mutate(priest_cnt = win_cnt / value)

                           return(priest_post)
                         },
                         otherwise = NA_character_)))

# make a plot
priest_p = comp_df %>%
  select(-dam_cnts) %>%
  unnest(dabom_est) %>%
  mutate(across(c(priest_cnt, priest_cnts),
                ~ . / 1000)) %>%
  ggplot(aes(x = priest_cnt,
             color = dam,
             fill = dam)) +
  geom_density(alpha = 0.2) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_vline(aes(xintercept = priest_cnts),
             lwd = 1,
             linetype = 2) +
  labs(x = "Estimate (x1000)",
       color = "Based on Counts at:",
       fill = "Based on Counts at:") +
  facet_wrap(~ year,
             scales = "free") +
  theme(legend.position = "bottom",
        axis.text.y = element_blank())

ggsave(here("analysis/figures",
            "PriestEstimates_AllYrs.pdf"),
       width = 8,
       height = 10)


# table with comparison of means
comp_diff = comp_df %>%
  select(-dam_cnts) %>%
  unnest(dabom_est) %>%
  group_by(year, dam, priest_cnts) %>%
  summarise(mean = mean(priest_cnt),
            median = median(priest_cnt),
            mode = estMode(priest_cnt),
            sd = sd(priest_cnt),
            skew = moments::skewness(priest_cnt),
            kurtosis = moments::kurtosis(priest_cnt),
            lowerCI = coda::HPDinterval(coda::as.mcmc(priest_cnt))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(priest_cnt))[,2],
            .groups = 'drop') %>%
  mutate(est_cv = sd / mean,
         in_ci = if_else(lowerCI <= priest_cnts &
                           upperCI >= priest_cnts,
                         T, F)) %>%
  mutate(dam = factor(dam,
                      levels = c("PriestRapids",
                                 "RockIsland",
                                 "RockyReach",
                                 "Wells",
                                 "Tumwater")))

comp_tab = comp_df %>%
  select(year,
         dam_cnts) %>%
  unnest(dam_cnts) %>%
  mutate(dam = factor(dam,
                      levels = c("PriestRapids",
                                 "Wanapum",
                                 "RockIsland",
                                 "RockyReach",
                                 "Wells",
                                 "Tumwater"))) %>%
  pivot_wider(names_from = dam,
              values_from = win_cnt) %>%
left_join(comp_diff %>%
    select(year, priest_cnts,
           dam,
           mean) %>%
    pivot_wider(names_from = dam,
                values_from = mean,
                names_prefix = "POM_")) %>%
  left_join(comp_diff %>%
    mutate(across(c(mean, median, mode),
                  ~ . - priest_cnts)) %>%
    select(year, priest_cnts,
           dam,
           mean) %>%
    pivot_wider(names_from = dam,
                values_from = mean,
                names_prefix = "diff_")) %>%
  left_join(comp_diff %>%
              mutate(across(c(mean, median, mode),
                            ~ (. - priest_cnts) / priest_cnts)) %>%
              select(year, priest_cnts,
                     dam,
                     mean) %>%
              pivot_wider(names_from = dam,
                          values_from = mean,
                          names_prefix = "perc_diff_")) %>%
  left_join(comp_diff %>%
              select(year, priest_cnts,
                     dam,
                     in_ci) %>%
              pivot_wider(names_from = dam,
                          values_from = in_ci,
                          names_prefix = "in_CI_")) %>%
  select(-priest_cnts)

write_csv(comp_tab,
          here("outgoing/other",
               "PriestCountComparison.csv"))

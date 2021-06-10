# Author: Kevin See
# Purpose: summarize DABOM results
# Created: 4/1/20
# Last Modified: 6/9/21
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(DABOM)
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(readxl)
library(STADEM)
library(writexl)
library(msm)
library(moments)
library(coda)
library(here)

#-----------------------------------------------------------------
# set year
yr = 2020

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data',
          'site_config.rda'))

# load compressed detections and biological data
load(here('analysis/data/derived_data/PITcleanr',
          paste0('UC_Steelhead_', yr, '.rda')))

# load JAGS MCMC results
load(here("analysis/data/derived_data/model_fits",
          paste0('PRA_DABOM_Steelhead_', yr,'.rda')))


# estimate final spawning location
tag_summ = summarizeTagData(filter_obs,
                            bio_df %>%
                              group_by(tag_code) %>%
                              slice(1) %>%
                              ungroup())

# look at which branch each tag was assigned to for spawning
brnch_df = buildNodeOrder(addParentChildNodes(parent_child, configuration)) %>%
  separate(col = path,
           into = paste("step", 1:max(.$node_order), sep = "_"),
           remove = F) %>%
  mutate(group = if_else(node == "PRA",
                         "Start",
                         if_else(grepl('LWE', path) | node %in% c("CLK"),
                                 "Wenatchee",
                                 if_else(grepl("ENL", path),
                                         "Entiat",
                                         if_else(grepl("LMR", path),
                                                 "Methow",
                                                 if_else(grepl("OKL", path) | node %in% c("FST"),
                                                         "Okanogan",
                                                         if_else(step_2 != "RIA" & !is.na(step_2),
                                                                 "BelowPriest",
                                                                 if_else(node == "WEA",
                                                                         "WellsPool",
                                                                         "Other")))))))) %>%
  select(-starts_with("step")) %>%
  mutate(group = factor(group,
                        levels = c("Wenatchee",
                                   "Entiat",
                                   "Methow",
                                   "Okanogan",
                                   "BelowPriest",
                                   "WellsPool",
                                   "Start",
                                   "Other")))

tag_summ %<>%
  left_join(brnch_df,
            by = c("spawn_node" = "node"))



# summarize detection probabilities
detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   filter_ch = filter_obs)

# which sites had detection probabilities fixed at 0% or 100%
detect_summ %>%
  filter(sd == 0)

# look at all the other sites
detect_summ %>%
  filter(sd > 0) %>%
  # arrange(desc(n_tags))
  arrange(desc(sd))

# compile all movement probabilities, and multiply them appropriately
trans_df = compileTransProbs_PRA(dabom_mod,
                                 parent_child)

# # summarize transition probabilities
# trans_summ = trans_df %>%
#   group_by(origin, param) %>%
#   summarise(mean = mean(value),
#             median = median(value),
#             mode = estMode(value),
#             sd = sd(value),
#             skew = moments::skewness(value),
#             kurtosis = moments::kurtosis(value),
#             lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
#             upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2],
#             .groups = "drop") %>%
#   mutate(across(c(mean, median, mode, sd, matches('CI$')),
#                 ~ if_else(. < 0, 0, .)))

#-----------------------------------------------------------------
# total escapement past Priest, by origin
start_date = paste0(yr-1, '0601')
end_date = paste0(yr, '0531')


# start with PIT-tag based reascension data
org_escape = queryPITtagData(damPIT = 'PRA',
                             spp = "Steelhead",
                             start_date = start_date,
                             end_date = end_date) %>%
  mutate(SpawnYear = yr) %>%
  mutate(across(TagIdAscentCount,
                tidyr::replace_na,
                0)) %>%
  mutate(ReAscent = ifelse(TagIdAscentCount > 1, T, F)) %>%
  group_by(Species, SpawnYear, Date) %>%
  summarise(tot_tags = n_distinct(TagId),
            reascent_tags = n_distinct(TagId[ReAscent]),
            .groups = "drop") %>%
  group_by(Species, SpawnYear) %>%
  summarise(across(matches('tags'),
                   sum,
                   na.rm = T),
            .groups = "drop") %>%
  mutate(reasc_rate = reascent_tags / tot_tags,
         reasc_rate_se = sqrt(reasc_rate * (1 - reasc_rate) / tot_tags)) %>%
  # add window counts
  bind_cols(getWindowCounts(dam = 'PRD',
                            spp = "Steelhead",
                            start_date = start_date,
                            end_date = end_date) %>%
              summarise_at(vars(win_cnt),
                           list(sum),
                           na.rm = T) %>%
              select(tot_win_cnt = win_cnt)) %>%
  mutate(adj_win_cnt = tot_win_cnt * (1 - reasc_rate),
         adj_win_cnt_se = tot_win_cnt * reasc_rate_se) %>%
  bind_cols(bio_df %>%
              group_by(origin) %>%
              summarise(n_tags = n_distinct(tag_code)) %>%
              mutate(prop = n_tags / sum(n_tags),
                     prop_se = sqrt((prop * (1 - prop)) / sum(n_tags)))) %>%
  rowwise() %>%
  mutate(tot_escp = adj_win_cnt * prop,
         tot_escp_se = msm::deltamethod(~ x1 * x2,
                                        mean = c(adj_win_cnt, prop),
                                        cov = diag(c(adj_win_cnt_se, prop_se)^2))) %>%
  select(Species, SpawnYear, origin, matches('escp'))

# translate movement estimates to escapement
escape_post = trans_df %>%
  left_join(org_escape %>%
              mutate(origin = recode(origin,
                                     "H" = 2,
                                     "W" = 1)) %>%
              group_by(origin) %>%
              summarise(tot_esc_samp = map2(tot_escp,
                                            tot_escp_se,
                                            .f = function(x, y) {
                                              tibble(tot_escp = rnorm(max(trans_df$iter),
                                                                      mean = x,
                                                                      sd = y)) %>%
                                                mutate(iter = 1:n())
                                            })) %>%
              unnest(cols = tot_esc_samp)) %>%
  mutate(escp = value * tot_escp)

escape_summ = escape_post %>%
  group_by(origin, location = param) %>%
  summarise(mean = mean(escp),
            median = median(escp),
            mode = estMode(escp),
            sd = sd(escp),
            skew = moments::skewness(escp),
            kurtosis = moments::kurtosis(escp),
            lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2],
            .groups = 'drop') %>%
  mutate(across(c(mean, median, mode, sd, matches('CI$')),
                ~ if_else(. < 0, 0, .))) %>%
  mutate(across(c(mean, median, mode, sd, skew, kurtosis, matches('CI$')),
                round,
                digits = 2)) %>%
  arrange(desc(origin), location) %>%
  tibble::add_column(species = "Steelhead",
                     spawn_year = yr,
                     .before = 0)



#-----------------------------------------------------------------
# create some summaries of biological information

bio_summ = tag_summ %>%
  group_by(group,
           origin,
           sex,
           age) %>%
  summarise(n_tags = n_distinct(tag_code),
            mean_FL = mean(fork_length, na.rm = T),
            .groups = "drop") %>%
  full_join(expand(tag_summ,
                   group,
                   origin,
                   nesting(sex, age))) %>%
  mutate(across(n_tags,
                replace_na,
                0)) %>%
  group_by(group, origin) %>%
  mutate(tot_tags = sum(n_tags)) %>%
  ungroup() %>%
  mutate(prop = n_tags / tot_tags,
         prop_se = sqrt((prop * (1 - prop)) / (tot_tags))) %>%
  arrange(group, origin, sex, age)

# population (or branch) level summary of escapement
pop_summ = escape_post %>%
  filter(param %in% c('LWE',
                      'ENL',
                      'LMR',
                      'OKL',
                      'ICH', 'JD1', 'JDA', 'PRH', 'PRO', 'PRV', 'RSH', 'TMF',
                      'WEA_bb')) %>%
  mutate(param = recode(param,
                        'LWE' = 'Wenatchee',
                        'ENL' = 'Entiat',
                        'LMR' = 'Methow',
                        'OKL' = 'Okanogan',
                        'ICH' = 'BelowPriest',
                        'JD1' = 'BelowPriest',
                        'JDA' = 'BelowPriest',
                        'PRH' = 'BelowPriest',
                        'PRO' = 'BelowPriest',
                        'PRV' = 'BelowPriest',
                        'RSH' = 'BelowPriest',
                        'TMF' = 'BelowPriest',
                        'WEA_bb' = "WellsPool")) %>%
  mutate(param = factor(param,
                        levels = levels(brnch_df$group))) %>%
  mutate(origin = recode(origin,
                         '1' = 'W',
                         '2' = 'H')) %>%
  group_by(chain, iter, origin, param) %>%
  summarize(across(escp,
                   sum),
            .groups = "drop") %>%
  group_by(origin, group = param) %>%
  summarise(mean = mean(escp),
            median = median(escp),
            mode = estMode(escp),
            sd = sd(escp),
            skew = moments::skewness(escp),
            kurtosis = moments::kurtosis(escp),
            lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2],
            .groups = 'drop') %>%
  mutate(across(c(mean, median, mode, sd, matches('CI$')),
                ~ if_else(. < 0, 0, .))) %>%
  mutate(across(c(mean, median, mode, sd, skew, kurtosis, matches('CI$')),
                round,
                digits = 2)) %>%
  arrange(desc(origin), group) %>%
  tibble::add_column(species = "Steelhead",
                     spawn_year = yr,
                     .before = 0)

#------------------------------------------------------------
# estimates of abundance by population, origin, sex and age
full_summ = pop_summ %>%
  # filter(group != "WellsPool") %>%
  select(species, spawn_year, group,
         origin, escape = mean, escp_se = sd) %>%
  left_join(bio_summ) %>%
  rowwise() %>%
  mutate(est = escape * prop,
         est_se = deltamethod(~ x1 * x2,
                              mean = c(escape, prop),
                              cov = diag(c(escp_se, prop_se)^2))) %>%
  ungroup() %>%
  mutate(across(est,
                ~ as.integer(round(est)))) %>%
  select(species:origin, sex, age, n_tags, mean_FL, starts_with("prop"), starts_with("est"))

#------------------------------------------------------------
# sex proportions by origin
sex_origin_summ = pop_summ %>%
  # filter(group != "WellsPool") %>%
  select(species, spawn_year, group,
         origin, escape = mean, escp_se = sd) %>%
  left_join(tag_summ %>%
              group_by(group, origin, sex) %>%
              summarise(n_tags = n_distinct(tag_code[!is.na(sex)]),
                        .groups = "drop") %>%
              filter(!is.na(sex)) %>%
              full_join(expand(tag_summ,
                               group,
                               origin, sex)) %>%
              mutate(across(n_tags,
                            replace_na,
                            0)) %>%
              group_by(group, origin) %>%
              mutate(total_sexed = sum(n_tags)) %>%
              mutate(prop = n_tags / total_sexed,
                     prop_se = sqrt((prop * (1 - prop)) / total_sexed)) %>%
              ungroup() %>%
              arrange(group, origin, sex)) %>%
  rowwise() %>%
  mutate(est = escape * prop,
         est_se = deltamethod(~ x1 * x2,
                              mean = c(escape, prop),
                              cov = diag(c(escp_se, prop_se)^2))) %>%
  ungroup() %>%
  select(species:origin, sex, n_tags, total_sexed, starts_with('prop'), starts_with("est")) %>%
  arrange(species, spawn_year, group, origin, sex)

# sex proportions overall
sex_summ = pop_summ %>%
  # filter(group != "WellsPool") %>%
  group_by(species, spawn_year, group) %>%
  summarize(escape = sum(mean),
            escp_se = sqrt(sum(sd^2)),
            .groups = "drop") %>%
  left_join(tag_summ %>%
              group_by(group, sex) %>%
              summarise(n_tags = n_distinct(tag_code[!is.na(sex)]),
                        .groups = "drop") %>%
              filter(!is.na(sex)) %>%
              full_join(expand(tag_summ,
                               group,
                               sex)) %>%
              mutate(across(n_tags,
                            replace_na,
                            0)) %>%
              group_by(group) %>%
              mutate(total_sexed = sum(n_tags)) %>%
              mutate(prop = n_tags / total_sexed,
                     prop_se = sqrt((prop * (1 - prop)) / total_sexed)) %>%
              ungroup() %>%
              arrange(group, sex)) %>%
  rowwise() %>%
  mutate(est = escape * prop,
         est_se = deltamethod(~ x1 * x2,
                              mean = c(escape, prop),
                              cov = diag(c(escp_se, prop_se)^2))) %>%
  ungroup() %>%
  select(species:group, sex, n_tags, total_sexed, starts_with('prop'), starts_with("est")) %>%
  arrange(species, spawn_year, group, sex)

# combine sex summaries
sex_all_summ = sex_origin_summ %>%
  bind_rows(sex_summ %>%
              mutate(origin = "All")) %>%
  mutate(across(origin,
                factor,
                levels = c("W", 'H', "All"))) %>%
  arrange(group, origin, sex) %>%
  mutate(across(est,
                round))

#------------------------------------------------------------
# age proportions by origin
age_origin_summ = pop_summ %>%
  # filter(group != "WellsPool") %>%
  select(species, spawn_year, group,
         origin, escape = mean, escp_se = sd) %>%
  left_join(tag_summ %>%
              group_by(group, origin, age) %>%
              summarise(n_tags = n_distinct(tag_code[!is.na(age)]),
                        .groups = "drop") %>%
              filter(!is.na(age)) %>%
              full_join(expand(tag_summ,
                               group,
                               origin, age)) %>%
              mutate(across(n_tags,
                            replace_na,
                            0)) %>%
              group_by(group, origin) %>%
              mutate(total_aged = sum(n_tags)) %>%
              mutate(prop = n_tags / total_aged,
                     prop_se = sqrt((prop * (1 - prop)) / total_aged)) %>%
              ungroup() %>%
              arrange(group, origin, age)) %>%
  rowwise() %>%
  mutate(est = escape * prop,
         est_se = deltamethod(~ x1 * x2,
                              mean = c(escape, prop),
                              cov = diag(c(escp_se, prop_se)^2))) %>%
  ungroup() %>%
  select(species:origin, age, n_tags, total_aged, starts_with('prop'), starts_with("est")) %>%
  arrange(species, spawn_year, group, origin, age)

# age proportions overall
age_summ = pop_summ %>%
  # filter(group != "WellsPool") %>%
  group_by(species, spawn_year, group) %>%
  summarize(escape = sum(mean),
            escp_se = sqrt(sum(sd^2)),
            .groups = "drop") %>%
  left_join(tag_summ %>%
              group_by(group, age) %>%
              summarise(n_tags = n_distinct(tag_code[!is.na(age)]),
                        .groups = "drop") %>%
              filter(!is.na(age)) %>%
              full_join(expand(tag_summ,
                               group,
                               age)) %>%
              mutate(across(n_tags,
                            replace_na,
                            0)) %>%
              group_by(group) %>%
              mutate(total_aged = sum(n_tags)) %>%
              mutate(prop = n_tags / total_aged,
                     prop_se = sqrt((prop * (1 - prop)) / total_aged)) %>%
              ungroup() %>%
              arrange(group, age)) %>%
  rowwise() %>%
  mutate(est = escape * prop,
         est_se = deltamethod(~ x1 * x2,
                              mean = c(escape, prop),
                              cov = diag(c(escp_se, prop_se)^2))) %>%
  ungroup() %>%
  select(species:group, age, n_tags, total_aged, starts_with('prop'), starts_with("est")) %>%
  arrange(species, spawn_year, group, age)

# combine age summaries
age_all_summ = age_origin_summ %>%
  bind_rows(age_summ %>%
              mutate(origin = "All")) %>%
  mutate(across(origin,
                factor,
                levels = c("W", 'H', "All"))) %>%
  arrange(group, origin, age) %>%
  mutate(across(est,
                round))

#------------------------------------------------------------
# origin proportion based on tags observed in each branch/population
org_summ = tag_summ %>%
  filter(!group %in% c("Start",
                       "Other")) %>%
  group_by(group, origin) %>%
  summarise(n_tags = n_distinct(tag_code),
            .groups = "drop") %>%
  pivot_wider(names_from = "origin",
              values_from = "n_tags",
              values_fill = 0) %>%
  mutate(prop_H = H / (H + W),
         prop_W = 1 - prop_H,
         prop_se = sqrt((prop_W * (1 - prop_W)) / (W + H)))

# # origin proportion based on DABOM estimates of abundance in each branch/population
# org_summ = pop_summ %>%
#   select(species:mean, sd) %>%
#   pivot_wider(names_from = "origin",
#               values_from = c("mean", "sd")) %>%
#   mutate(prop_W = mean_W / (mean_W + mean_H),
#          prop_H = 1 - prop_W) %>%
#   rowwise() %>%
#   mutate(prop_se = msm::deltamethod(~ x1 / (x1 + x2),
#                                     mean = c(mean_W, mean_H),
#                                     cov = diag(c(sd_W, sd_H)^2))) %>%
#   arrange(group)

# biological summaries, based on tags detected within each stream / population
bio_list = list('Origin' = org_summ,
                'Sex' = sex_all_summ,
                "Sexed Tags" = sex_all_summ %>%
                  select(-starts_with(c("prop", "est"))) %>%
                  pivot_wider(names_from = "sex",
                              values_from = "n_tags",
                              values_fill = 0,
                              names_sort = T),
                'Age' = age_all_summ,
                'Aged Tags' = age_all_summ %>%
                  select(-starts_with(c("prop", "est"))) %>%
                  pivot_wider(names_from = "age",
                              values_from = "n_tags",
                              values_fill = 0,
                              names_sort = T),
                'BioSummary' = full_summ)


#-----------------------------------------------------------------
# write results to an Excel file
save_list = c(list('Population Escapement' = pop_summ %>%
                     select(-skew, -kurtosis) %>%
                     mutate_at(vars(mean:upperCI),
                               list(round),
                               digits = 1) %>%
                     rename(estimate = mean,
                            se = sd) %>%
                     select(-median, -mode),
                   'All Escapement' = escape_summ %>%
                     select(-skew, -kurtosis) %>%
                     mutate_at(vars(mean:upperCI),
                               list(round),
                               digits = 1) %>%
                     rename(estimate = mean,
                            se = sd) %>%
                     select(-median, -mode),
                   'Detection' = detect_summ %>%
                     mutate(across(-c(node, n_tags),
                                   round,
                                   digits = 3)) %>%
                     rename(estimate = mean,
                            se = sd) %>%
                     select(-median, -mode),
                   'Tag Summary' = tag_summ),
              bio_list)

writexl::write_xlsx(x = save_list,
                    path = here('outgoing/estimates',
                                paste0('UC_Steelhead_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx')))

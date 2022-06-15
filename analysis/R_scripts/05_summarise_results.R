# Author: Kevin See
# Purpose: summarize DABOM results
# Created: 4/1/20
# Last Modified: 6/15/2022
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
# load configuration and site_df data
load(here('analysis/data/derived_data',
          'site_config.rda'))

#-----------------------------------------------------------------
# set year
yr = 2021

#-----------------------------------------------------------------
# run for set of years
# for(yr in 2011:2021) {

  cat(paste("Working on", yr, "\n\n"))

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
                                   parent_child) %>%
    mutate(origin = recode(origin,
                           "2" = "H",
                           "1" = "W"))

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
    filter(!str_detect(TagId, "000.0")) %>%
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
    select(Species, SpawnYear, origin, reasc_rate, matches('escp'))


  # #----------------------------------------------------------------
  # # better way to use upstream dam counts
  # # add origin to prepped capture histories
  # pit_obs = prepped_ch %>%
  #   left_join(bio_df %>%
  #               select(tag_code,
  #                      origin)) %>%
  #   select(tag_code, origin,
  #          everything())
  #
  # # get the total dam counts from various dams
  # dam_cnts = tibble(dam = c("PriestRapids",
  #                           "RockIsland",
  #                           "RockyReach",
  #                           "Wells",
  #                           "Tumwater"),
  #                   dam_code = c("PRD",
  #                                "RIS",
  #                                "RRH",
  #                                "WEL",
  #                                "TUM"),
  #                   year = yr) %>%
  #   rowwise() %>%
  #   mutate(win_cnt = map_dbl(dam_code,
  #                             .f = function(x, y) {
  #                               STADEM::getWindowCounts(dam = x,
  #                                                       spp = "Steelhead",
  #                                                       start_date = start_date,
  #                                                       end_date = end_date) %>%
  #                                 summarise_at(vars(win_cnt),
  #                                              list(sum),
  #                                              na.rm = T) %>%
  #                                 pull(win_cnt)
  #                             })) %>%
  #   ungroup() %>%
  #   select(year, everything())
  #
  # # compile and calculate Priest equivalents
  # pit_move_df = tibble(dam = c("PriestRapids",
  #                              "RockIsland",
  #                              "RockyReach",
  #                              "Wells"),
  #                      dam_code = c("PRD",
  #                                   "RIS",
  #                                   "RRH",
  #                                   "WEL"),
  #                      pit_code = c("PRA",
  #                                   "RIA",
  #                                   "RRF",
  #                                   "WEA")) %>%
  #   add_column(year = yr,
  #              .before = 0) %>%
  #   crossing(origin = c("W", "H")) %>%
  #   left_join(pit_obs %>%
  #               group_by(origin) %>%
  #               summarize(n_tags = n_distinct(tag_code))) %>%
  #   mutate(n_obs = map2_int(pit_code,
  #                           origin,
  #                           .f = function(x, y) {
  #                             pit_obs %>%
  #                               filter(origin == y) %>%
  #                               summarize(n_obs = n_distinct(tag_code[node == x])) %>%
  #                               pull(n_obs)
  #                           })) %>%
  #   mutate(n_upstrm = map2_int(pit_code,
  #                              origin,
  #                              .f = function(x, y) {
  #                                pit_obs %>%
  #                                  filter(origin == y) %>%
  #                                  summarize(n_path = n_distinct(tag_code[str_detect(path, paste0(" ", x))])) %>%
  #                                  pull(n_path)
  #                              }),
  #          n_upstrm = if_else(dam == "PriestRapids",
  #                             n_obs,
  #                             n_upstrm)) %>%
  #   # generate proportion hatchery/wild based on upstream tags
  #   group_by(dam) %>%
  #   mutate(prop_org = n_upstrm / sum(n_upstrm),
  #          prop_org_se = sqrt((prop_org * (1 - prop_org)) / sum(n_upstrm))) %>%
  #   # generate proportion hatchery/wild based on tags detected at dam
  #   # mutate(prop_org = n_obs / sum(n_obs),
  #   #        prop_org_se = sqrt((prop_org * (1 - prop_org)) / sum(n_obs))) %>%
  #   ungroup() %>%
  #   left_join(detect_summ %>%
  #               select(node, mean, sd),
  #             by = c("pit_code" = "node")) %>%
  #   mutate(est_tags = n_obs / mean,
  #          est_tags_se = n_obs * sd / mean^2,
  #          trans_est = est_tags / n_tags,
  #          trans_se = est_tags_se / n_tags) %>%
  #   mutate(trans_est = if_else(dam == "PriestRapids",
  #                              1,
  #                              trans_est),
  #          trans_se = if_else(dam == "PriestRapids",
  #                             0,
  #                             trans_se)) %>%
  #   left_join(dam_cnts) %>%
  #   rowwise() %>%
  #   mutate(priest_cnt = win_cnt * prop_org / trans_est,
  #          priest_cnt_se = msm::deltamethod(~ x1 * x2 / x3,
  #                                           mean = c(win_cnt,
  #                                                    prop_org,
  #                                                    trans_est),
  #                                           cov = diag(c(0,
  #                                                        prop_org_se,
  #                                                        trans_se)^2))) %>%
  #   ungroup() %>%
  #   # add re-ascension rate at Priest
  #   left_join(queryPITtagData(damPIT = 'PRA',
  #                             spp = "Steelhead",
  #                             start_date = start_date,
  #                             end_date = end_date) %>%
  #               filter(!str_detect(TagId, "000.0")) %>%
  #               mutate(SpawnYear = yr) %>%
  #               mutate(across(TagIdAscentCount,
  #                             tidyr::replace_na,
  #                             0)) %>%
  #               mutate(ReAscent = ifelse(TagIdAscentCount > 1, T, F)) %>%
  #               mutate(origin = fct_recode(RearType,
  #                                          "W" = "U")) %>%
  #               group_by(Species, SpawnYear, Date, origin) %>%
  #               summarise(tot_tags = n_distinct(TagId),
  #                         reascent_tags = n_distinct(TagId[ReAscent]),
  #                         .groups = "drop") %>%
  #               group_by(Species, SpawnYear, origin) %>%
  #               summarise(across(matches('tags'),
  #                                sum,
  #                                na.rm = T),
  #                         .groups = "drop") %>%
  #               mutate(reasc_rate = reascent_tags / tot_tags,
  #                      reasc_rate_se = sqrt(reasc_rate * (1 - reasc_rate) / tot_tags)) %>%
  #               select(origin, starts_with("reasc_rate"))) %>%
  #   rowwise() %>%
  #   mutate(tot_escp = priest_cnt * (1 - reasc_rate),
  #          tot_escp_se = msm::deltamethod(~ x1 * (1 - x2),
  #                                         mean = c(priest_cnt,
  #                                                  reasc_rate),
  #                                         cov = diag(c(priest_cnt_se,
  #                                                      reasc_rate_se)^2))) %>%
  #   ungroup()
  #
  # org_escape <- pit_move_df %>%
  #   filter(dam == "PriestRapids") %>%
  #   # filter(dam == "RockIsland") %>%
  #   mutate(Species = "Steelhead") %>%
  #   select(Species, SpawnYear = year,
  #          origin,
  #          starts_with("priest_cnt"),
  #          reasc_rate, matches('escp'))


  # translate movement estimates to escapement
  escape_post = trans_df %>%
    left_join(org_escape %>%
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
              sd_FL = sd(fork_length, na.rm = T),
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
                  'Biological Summary' = full_summ)


  #-----------------------------------------------------------------
  # break down hatchery estimates by mark type
  #-----------------------------------------------------------------
  mark_tag_summ = tag_summ %>%
    # filter(origin == "H") %>%
    mutate(cwt = if_else(is.na(cwt), F,
                         if_else(cwt %in% c("SN", "BD"),
                                 T, NA)),
           ad_clip = if_else(is.na(ad_clip), F,
                             if_else(ad_clip == "AD", T, NA))) %>%
    mutate(ad_clip_chr = recode(as.character(ad_clip),
                                "TRUE" = "AD",
                                "FALSE" = "AI"),
           cwt_chr = recode(as.character(cwt),
                            "TRUE" = "CWT",
                            "FALSE" = "noCWT")) %>%
    tidyr::unite("mark_grp", ad_clip_chr, cwt_chr, remove = T) %>%
    mutate(mark_grp = if_else(origin == "W",
                              "Wild",
                              mark_grp))

  # proportions of each type of fish (H/W, Ad-clip, CWT combinations) past each site
  mark_grp_prop = mark_tag_summ %>%
    mutate(spawn_site = str_remove(spawn_node, "B0$"),
           spawn_site = str_remove(spawn_site, "A0$"),
           spawn_site = recode(spawn_site,
                               "S" = "SA0")) %>%
    left_join(buildPaths(parent_child) %>%
                separate(path,
                         into = paste("site", 1:8, sep = "_")) %>%
                pivot_longer(starts_with('site_'),
                             names_to = "node_order",
                             values_to = "site_code") %>%
                filter(!is.na(site_code)) %>%
                select(-node_order),
              by = c("spawn_site" = "end_loc")) %>%
    group_by(origin, site_code, ad_clip, cwt, mark_grp) %>%
    summarise(n_tags = n_distinct(tag_code),
              .groups = "drop") %>%
    right_join(crossing(site_code = union(parent_child$parent, parent_child$child),
                        expand(mark_tag_summ, nesting(origin, ad_clip, cwt, mark_grp)))) %>%
    arrange(site_code, mark_grp) %>%
    mutate(across(n_tags,
                  replace_na,
                  0)) %>%
    group_by(origin, site_code) %>%
    mutate(tot_tags = sum(n_tags),
           prop = n_tags / tot_tags,
           # using normal approximation
           prop_se = sqrt((prop * (1 - prop))/tot_tags)) %>%
    ungroup()

  # generate posterior samples of mark proportions
  set.seed(6)
  n_iter = max(escape_post$iter)
  prop_samps = mark_grp_prop %>%
    filter(tot_tags > 0) %>%
    group_by(origin, site_code,
             tot_tags) %>%
    nest() %>%
    mutate(samp = map(data,
                      .f = function(x) {
                        rmultinom(n_iter, sum(x$n_tags), x$prop) %>%
                          set_colnames(1:ncol(.)) %>%
                          set_rownames(x$mark_grp) %>%
                          as_tibble(rownames = 'mark_grp') %>%
                          pivot_longer(-1,
                                       names_to = "iter",
                                       values_to = "n_tags") %>%
                          mutate(across(iter,
                                        as.integer)) %>%
                          group_by(iter) %>%
                          mutate(prop = n_tags / sum(n_tags)) %>%
                          arrange(iter, mark_grp) %>%
                          ungroup()
                      })) %>%
    ungroup() %>%
    select(-data) %>%
    unnest(samp) %>%
    left_join(mark_grp_prop %>%
                select(-n_tags,
                       -starts_with("prop"))) %>%
    bind_rows(mark_grp_prop %>%
                filter(tot_tags == 0) %>%
                crossing(iter = 1:n_iter) %>%
                select(-prop_se)) %>%
    arrange(site_code, origin, mark_grp, iter) %>%
    select(iter, any_of(names(mark_grp_prop)))

  # posterior samples
  mark_post = escape_post %>%
    inner_join(prop_samps,
               by = c("iter", "origin",
                      "param" = "site_code")) %>%
    mutate(across(c(prop,
                    n_tags,
                    tot_tags),
                  replace_na,
                  0)) %>%
    mutate(prop = if_else(value == 0,
                          0,
                          prop)) %>%
    mutate(n_fish = escp * prop)

  mark_grp_summ = mark_post %>%
    group_by(origin,
             location = param,
             ad_clip,
             cwt,
             mark_grp) %>%
    summarise(across(n_fish,
                     list(mean = mean,
                          median = median,
                          mode = estMode,
                          se = sd,
                          skew = moments::skewness,
                          kurtosis = moments::kurtosis,
                          lowerCI = ~ coda::HPDinterval(coda::as.mcmc(.x))[,1],
                          upperCI = ~ coda::HPDinterval(coda::as.mcmc(.x))[,2]),
                     na.rm = T,
                     .names = "{.fn}"),
              .groups = "drop") %>%
    mutate(across(mode,
                  ~ if_else(. < 0, 0, .))) %>%
    left_join(mark_grp_prop %>%
                rename(location = site_code,
                       proportion = prop)) %>%
    select(origin:mark_grp,
           n_tags:proportion,
           prop_se,
           everything()) %>%
    arrange(location, mark_grp)

  pop_mark_grp = mark_grp_summ %>%
    filter(location %in% c("LWE", 'ENL', 'LMR', 'OKL')) %>%
    mutate(population = recode(location,
                               'LWE' = 'Wenatchee',
                               'ENL' = 'Entiat',
                               'LMR' = 'Methow',
                               'OKL' = 'Okanogan'),
           population = factor(population,
                               levels = c('Wenatchee',
                                          'Entiat',
                                          'Methow',
                                          'Okanogan'))) %>%
    arrange(population, mark_grp) %>%
    select(origin, population, everything(), -location)

  # biological summary that includes mark group
  mark_grp_bio_summ = mark_tag_summ %>%
    group_by(group,
             mark_grp,
             origin,
             sex,
             age) %>%
    summarise(n_tags = n_distinct(tag_code),
              mean_FL = mean(fork_length, na.rm = T),
              sd_FL = sd(fork_length, na.rm = T),
              .groups = "drop") %>%
    full_join(expand(mark_tag_summ,
                     group,
                     nesting(mark_grp, origin),
                     nesting(sex, age))) %>%
    mutate(across(n_tags,
                  replace_na,
                  0)) %>%
    group_by(group, origin, mark_grp) %>%
    mutate(tot_tags = sum(n_tags)) %>%
    ungroup() %>%
    mutate(prop = n_tags / tot_tags,
           prop_se = sqrt((prop * (1 - prop)) / (tot_tags))) %>%
    arrange(group, origin, mark_grp, sex, age)

  mark_grp_list = list('Mark Group Population' = pop_mark_grp %>%
                         select(-skew, -kurtosis) %>%
                         mutate_at(vars(mean:mode),
                                   list(round),
                                   digits = 0) %>%
                         mutate_at(vars(se:upperCI),
                                   list(round),
                                   digits = 1) %>%
                         rename(estimate = mean) %>%
                         select(-median, -mode),
                       'Mark Group Site Escapement' = mark_grp_summ %>%
                         select(-skew, -kurtosis) %>%
                         mutate_at(vars(mean:mode),
                                   list(round),
                                   digits = 0) %>%
                         mutate_at(vars(se:upperCI),
                                   list(round),
                                   digits = 1) %>%
                         rename(estimate = mean) %>%
                         select(-median, -mode),
                       'Mark Group Tag Summary' = mark_tag_summ,
                       'Mark Group Bio Summary' = mark_grp_bio_summ)

  #-----------------------------------------------------------------
  # write results to an Excel file
  save_list = c(list('Population Escapement' = pop_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate(across(mean:mode,
                                     janitor::round_half_up)) %>%
                       mutate(across(sd:upperCI,
                                     round,
                                     digits = 1)) %>%
                       rename(estimate = mean,
                              se = sd) %>%
                       select(-median, -mode),
                     'All Escapement' = escape_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate(across(mean:mode,
                                     janitor::round_half_up)) %>%
                       mutate(across(sd:upperCI,
                                     round,
                                     digits = 1)) %>%
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
                bio_list,
                mark_grp_list)

  writexl::write_xlsx(x = save_list,
                      path = here('outgoing/estimates',
                                  paste0('UC_Steelhead_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx')))

# }



#---------------------------------------------------
# compare estimates with some dam counts
#---------------------------------------------------

# query dam counts at Priest, Rock Island, Rocky Reach, Wells and Tumwater
dam_cnts = tibble(dam = c("RockIsland",
                          "RockyReach",
                          "Wells",
                          "Tumwater"),
                  dam_code = c("RIS",
                               "RRH",
                               "WEL",
                               "TUM")) %>%
  crossing(spp = c("Steelhead",
                   "Wild_Steelhead")) %>%
  mutate(win_cnt = map2_dbl(dam_code,
                            spp,
                            .f = function(x, y) {
                              STADEM::getWindowCounts(dam = x,
                                                      spp = y,
                                                      start_date = start_date,
                                                      end_date = end_date) %>%
                                summarise_at(vars(win_cnt),
                                             list(sum),
                                             na.rm = T) %>%
                                pull(win_cnt)
                            })) %>%
  pivot_wider(names_from = spp,
              values_from = win_cnt) %>%
  rename(All = Steelhead,
         W = Wild_Steelhead) %>%
  mutate(H = All - W) %>%
  pivot_longer(cols = c("All", "W", 'H'),
               names_to = 'origin',
               values_to = "win_cnt")

dam_est = escape_summ %>%
  filter(location %in% c("RIA",
                         "RRF",
                         "WEA",
                         "TUM")) %>%
  mutate(dam = recode(location,
                      "RIA" = "RockIsland",
                      "RRF" = "RockyReach",
                      "WEA" = "Wells",
                      "TUM" = "Tumwater"))
comp_df = dam_est %>%
  group_by(dam) %>%
  summarize(across(c(mean:mode, ends_with("CI")),
                   sum),
            .groups = "drop") %>%
  mutate(origin = "All") %>%
  bind_rows(dam_est %>%
              select(dam,
                     origin,
                     mean:mode,
                     ends_with("CI"))) %>%
  left_join(dam_cnts) %>%
  rowwise() %>%
  mutate(larger = if_else(win_cnt > mean,
                          "Counts",
                          "DABOM")) %>%
  mutate(in_ci = between(win_cnt, lowerCI, upperCI)) %>%
  ungroup() %>%
  mutate(dam = fct_relevel(dam,
                           "Tumwater",
                           after = Inf)) %>%
  arrange(dam,
          origin)

# add adjustments if re-ascension rate is the same everywhere as it is at Priest
comp_df %<>%
  add_column(reasc_rate = unique(org_escape$reasc_rate)) %>%
  mutate(adj_win_cnt = win_cnt * (1 - reasc_rate))

comp_df %>%
  ggplot(aes(x = dam,
             y = mean,
             color = "DABOM")) +
  geom_errorbar(aes(ymin = lowerCI,
                    ymax = upperCI),
                width = 0) +
  geom_point(size = 3) +
  geom_point(aes(y = win_cnt,
                 color = "Dam Count"),
             size = 3,
             position = position_dodge(width = 1)) +
  geom_point(aes(y = adj_win_cnt,
                 color = "Adj. Dam Count"),
             size = 3,
             position = position_dodge(width = 1)) +
  scale_color_manual(values = c("DABOM" = "gray20",
                                "Dam Count" = "red",
                                "Adj. Dam Count" = "blue"),
                     name = "Source") +
  facet_wrap(~ origin,
             scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  labs(x = "Dam",
       y = "Estimate",
       title = paste("Steelhead", yr),
       subtitle = "Using Rock Island Counts")

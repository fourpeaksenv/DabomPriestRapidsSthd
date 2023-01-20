# Author: Kevin See
# Purpose: format results to be saved
# Created: 11/30/22
# Last Modified: 12/21/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(UCSthdReddObsErr)
library(janitor)
library(writexl)
library(magrittr)
library(msm)
library(here)

# load age table
data("age_table")

#-----------------------------------------------------------------
# function to fix all table names
makeTableNms = function(df) {
  df |>
    rlang::set_names(nm = function(x) {
      x |>
        str_replace_all("_", " ") |>
        str_to_title() |>
        str_replace("Nos", "NOS") |>
        str_replace("Hos", "HOS") |>
        str_replace("Se$", "SE") |>
        str_replace("Sd$", "SD") |>
        str_replace("Cv$", "CV") |>
        str_replace("Fpr", "FpR") |>
        str_replace("Phos", "pHOS") |>
        str_replace("Lowerci", "LCI") |>
        str_replace("Upperci", "UCI") |>
        str_replace("Cwt$", "CWT")
    }) %>%
    return()
}


#-----------------------------------------------------------------
# generate tables of redds and spawners

# set some thresholds
# minimum number of total redds observed
min_redds = 2
# minimum number of weeks with at least one new redd observed
min_non0_wks = 3

# redd_spwn = NULL
# spwn_df = NULL

spwn_est <- crossing(population = c("Wenatchee",
                                    "Methow"),
                     spawn_year = c(2014:2022)) %>%
  mutate(run_year = spawn_year -1) %>%
  select(run_year, spawn_year, population) %>%
  mutate(file_nm = map2_chr(population,
                            spawn_year,
                            .f = function(pop_nm,
                                          yr) {
                              file_nm = here('O:/Documents/Git/MyProjects/UCSthdReddObsErr/analysis/data/derived_data',
                                             paste0(if_else(pop_nm == "Wenatchee",
                                                            'wen_',
                                                            if_else(pop_nm == "Methow",
                                                                    "met_",
                                                                    NA_character_)), yr, '.rda'))
                              return(file_nm)
                            }),
         file_exists = map_lgl(file_nm,
                               .f = file.exists)) %>%
  filter(file_exists) %>%
  mutate(results_list = map2(file_nm,
                             population,
                            .f = possibly(function(file_name,
                                          pop_nm) {
                              load(file_name)

                              # for 2020, we'll need to do something different
                              if(str_detect(file_name, "_2020.rda",
                                            negate = T)) {

                              results_lst <- redd_df |>
                                bind_rows(redds_below_arrays) |>
                                group_by(Index, SurveyType) |>
                                nest(redd_data = -c(Index, SurveyType)) |>
                                filter(!(is.na(Index) & is.na(SurveyType))) |>
                                mutate(summ_lst = map2(redd_data, Index,
                                                       .f = function(x, y) {
                                                         summarizeRedds(redd_df = x,
                                                                        group_vars = c("Reach"),
                                                                        summ_vars = c("River",
                                                                                      "Location"),
                                                                        min_non0_wks = min_non0_wks,
                                                                        min_redds = min_redds,
                                                                        use_cor = if_else(y == "Y", T, F))
                                                       })) |>
                                ungroup()

                              redd_results <- results_lst |>
                                mutate(rch_res = map(summ_lst,
                                                     "rch_est")) |>
                                select(Index,
                                       SurveyType,
                                       rch_res) |>
                                unnest(rch_res) |>
                                ungroup() |>
                                mutate(Reach = factor(Reach,
                                                      levels = c(paste0('W', 1:10),
                                                                 "C1",
                                                                 "N1",
                                                                 "P1",
                                                                 paste0("MRW", 1:8)))) |>
                                arrange(Reach,
                                        Index) |>
                                relocate(River, Reach, Location,
                                         .before = 1)

                              spwn_rch = redd_results |>
                                select(River, Location, Reach, Index, redd_est, redd_se) |>
                                left_join(fpr_df |>
                                            select(Location, starts_with('fpr')),
                                          by = "Location") |>
                                left_join(fpr_df |>
                                            select(Location, starts_with('phos')),
                                          by = "Location") |>
                                rowwise() |>
                                mutate(Tot_Spawners = redd_est * fpr,
                                       Hatchery = redd_est * fpr * phos,
                                       Natural = redd_est * fpr * (1 - phos),
                                       Hatchery_SE = ifelse(Tot_Spawners > 0,
                                                            deltamethod(~ x1 * x2 * x3,
                                                                        mean = c(redd_est, fpr, phos),
                                                                        cov = diag(c(redd_se, fpr_se, phos_se)^2)),
                                                            NA),
                                       Natural_SE = ifelse(Tot_Spawners > 0,
                                                           deltamethod(~ x1 * x2 * (1 - x3),
                                                                       mean = c(redd_est, fpr, phos),
                                                                       cov = diag(c(redd_se, fpr_se, phos_se)^2)),
                                                           NA)) |>
                                ungroup() |>
                                select(River, Location, Reach,
                                       Type = Index,
                                       Hatchery:Natural_SE) |>
                                mutate(Type = recode(Type,
                                                     'Y' = 'Index',
                                                     'N' = 'Non-Index')) |>
                                bind_rows(trib_spawners |>
                                            pivot_wider(names_from = Origin,
                                                        values_from = c(Spawners,
                                                                        Spawners_SE),
                                                        names_glue = "{Origin}_{.value}") |>
                                            rlang::set_names(function(x) str_remove(x, "_Spawners")) |>
                                            mutate(Location = recode(Location,
                                                                     'CHL' = 'Chiwawa',
                                                                     'CHM' = 'Chumstick',
                                                                     'CHW' = 'Chiwaukum',
                                                                     'ICL' = 'Icicle',
                                                                     'LWN' = 'Little Wenatchee',
                                                                     'MCL' = 'Mission',
                                                                     'NAL' = 'Nason',
                                                                     'PES' = 'Peshastin',
                                                                     'WTL' = 'White River',
                                                                     "GLC" = "Gold",
                                                                     "LBC" = "Libby",
                                                                     "MSH" = "Methow Fish Hatchery",
                                                                     "MRW" = "Upper Methow",
                                                                     "TWR" = "Twisp",
                                                                     "CRW" = "Chewuch",
                                                                     "SCP" = "Spring Creek",
                                                                     "BVC" = "Beaver")) |>
                                            mutate(River = Location) |>
                                            mutate(Type = 'DABOM')) |>
                                mutate(Location = factor(Location, levels = c('Below_TUM',
                                                                              'TUM_bb',
                                                                              'Icicle',
                                                                              'Peshastin',
                                                                              'Mission',
                                                                              'Chumstick',
                                                                              'Chiwaukum',
                                                                              'Chiwawa',
                                                                              'Nason',
                                                                              'Little Wenatchee',
                                                                              'White River',
                                                                              "Lower Methow",
                                                                              "Gold",
                                                                              "Libby",
                                                                              "Methow Fish Hatchery",
                                                                              "Upper Methow",
                                                                              "Twisp",
                                                                              "Chewuch",
                                                                              "Spring Creek",
                                                                              "Beaver"))) |>
                                arrange(Location, Reach, Type) |>
                                mutate(h_cv = Hatchery_SE / Hatchery,
                                       n_cv = Natural_SE / Natural)


                              redd_spwn_yr <- redd_results |>
                                mutate(reach_type = recode(Index,
                                                           'Y' = 'Index',
                                                           'N' = 'Non-Index')) |>
                                left_join(fpr_df %>%
                                            select(Location,
                                                   starts_with("fpr")),
                                          by = "Location") |>
                                select(stream = River,
                                       Reach,
                                       reach_type,
                                       observed_redds = tot_feat,
                                       redd_est,
                                       redd_se,
                                       fpr,
                                       fpr_se) %>%
                                mutate(redd_cv = redd_se / redd_est) %>%
                                relocate(redd_cv,
                                         .before = fpr) %>%
                                left_join(spwn_rch %>%
                                            select(stream = River,
                                                   Reach,
                                                   reach_type = Type,
                                                   NOS_est = Natural,
                                                   NOS_se = Natural_SE,
                                                   NOS_cv = n_cv,
                                                   HOS_est = Hatchery,
                                                   HOS_se = Hatchery_SE,
                                                   HOS_cv = h_cv),
                                          by = c("stream", "Reach", "reach_type"))

                              spwn_yr <- spwn_rch |>
                                group_by(River) |>
                                summarize(across(c(nos_est = Natural,
                                                   hos_est = Hatchery),
                                                 sum,
                                                 na.rm = T),
                                          across(c(nos_se = Natural_SE,
                                                   hos_se = Hatchery_SE),
                                                 ~ sqrt(sum(.^2, na.rm = T))),
                                          .groups = "drop") %>%
                                mutate(nos_cv = nos_se / nos_est,
                                       hos_cv = hos_se / hos_est) %>%
                                mutate(total_spawners = nos_est + hos_est) |>
                                mutate(across(River,
                                              recode,
                                              "Wenatchee" = "Wenatchee (mainstem)",
                                              "Methow" = "Methow (mainstem)"),
                                       across(River,
                                              as.factor),
                                       across(River,
                                              fct_relevel,
                                              "Wenatchee (mainstem)",
                                              "Methow (mainstem)",
                                              after = Inf)) |>
                                arrange(River) |>
                                select(stream = River,
                                       total_spawners,
                                       starts_with("nos"),
                                       starts_with("hos"))

                              spwn_yr %<>%
                                bind_rows(spwn_yr %>%
                                          summarize(stream = "Total",
                                                    across(c(total_spawners,
                                                             nos_est,
                                                             hos_est),
                                                           sum),
                                                    across(ends_with("_se"),
                                                           ~ sqrt(sum(.^2, na.rm = T))),
                                                    nos_cv = nos_se / nos_est,
                                                    hos_cv = hos_se / hos_est,
                                                    .groups = "drop")) %>%
                                rowwise() %>%
                                mutate(phos = hos_est / (hos_est + nos_est),
                                       phos_se = deltamethod(~ x1 / (x1 + x2),
                                                             mean = c(hos_est,
                                                                      nos_est),
                                                             cov = diag(c(hos_se,
                                                                          nos_se)^2))) %>%
                                ungroup() %>%
                                mutate(across(stream,
                                              as.factor),
                                       across(stream,
                                              fct_relevel,
                                              c("Wenatchee (mainstem)",
                                                "Methow (mainstem)",
                                                "Total"),
                                              after = Inf)) %>%
                                arrange(stream)


                              list(results_lst = results_lst,
                                   fpr = fpr_df,
                                   redd_results = redd_results,
                                   spwn_rch = spwn_rch,
                                   redd_spwn_yr = redd_spwn_yr,
                                   spwn_yr = spwn_yr) %>%
                                return()

                              } else {
                                # data from WDFW (Nate Fuchs' radio telemetry study)
                                rt_df = tibble(year = rep(2015:2016, each = 2),
                                               Origin = rep(c("Hatchery", "Natural"), 2),
                                               ow_fish = c(20, 25, 4, 12),
                                               surv_fish = c(16, 24, 3, 9)) %>%
                                  mutate(phi = surv_fish / ow_fish,
                                         phi_se = sqrt((phi * (1 - phi))/ ow_fish))

                                # add years together
                                rt_df %<>%
                                  bind_rows(rt_df %>%
                                              group_by(Origin) %>%
                                              summarize(across(c(ow_fish, surv_fish),
                                                               sum),
                                                        .groups = "drop") %>%
                                              mutate(phi = surv_fish / ow_fish,
                                                     phi_se = sqrt((phi * (1 - phi))/ ow_fish))) %>%
                                  mutate(across(year,
                                                as.factor)) %>%
                                  mutate(across(year,
                                                fct_explicit_na,
                                                na_level = "Total"))

                                data("removal_df")
                                yr = str_extract(file_name, "[:digit:]+")

                                rem_df = removal_df %>%
                                  filter(Year == yr,
                                         Subbasin == pop_nm) %>%
                                  select(-Year)

                                spwn_yr <- escp_wen %>%
                                  filter(Area == "Wen_all") %>%
                                  left_join(rem_df %>%
                                              group_by(Origin) %>%
                                              summarize(across(rem,
                                                               sum),
                                                        .groups = "drop"),
                                            by = "Origin") %>%
                                  rowwise() %>%
                                  mutate(escp = max(0, estimate - rem)) %>%
                                  ungroup() %>%
                                  full_join(rt_df %>%
                                              group_by(Origin) %>%
                                              summarize(across(c(ow_fish, surv_fish),
                                                               sum),
                                                        .groups = "drop") %>%
                                              mutate(phi = surv_fish / ow_fish,
                                                     phi_se = sqrt((phi * (1 - phi))/ ow_fish)) %>%
                                              select(-ends_with("fish")),
                                            by = "Origin") %>%
                                  rowwise() %>%
                                  mutate(all_spwn = escp * phi) %>%
                                  mutate(all_spwn_se = msm::deltamethod(~ x1 * x2,
                                                                        mean = c(escp, phi),
                                                                        cov = diag(c(se, phi_se)^2))) %>%
                                  ungroup() %>%
                                  left_join(trib_spawners %>%
                                              group_by(Origin) %>%
                                              summarize(across(Spawners,
                                                               sum),
                                                        across(Spawners_SE,
                                                               ~ sqrt(sum(.^2))),
                                                        .groups = "drop") %>%
                                              rename(trib_spwn = Spawners,
                                                     trib_se = Spawners_SE),
                                            by = "Origin") %>%
                                  mutate(main_spwn = all_spwn - trib_spwn,
                                         main_spwn_se = sqrt(all_spwn_se^2 + trib_se^2),
                                         across(main_spwn,
                                                ~ if_else(. < 0, 0, .))) %>%
                                  mutate(Location = "Wenatchee (mainstem)",
                                         River = Location,
                                         Type = "RT") %>%
                                  select(River,
                                         Origin,
                                         Location,
                                         Type,
                                         Spawners = main_spwn,
                                         Spawners_SE = main_spwn_se) %>%
                                  bind_rows(trib_spawners %>%
                                              mutate(Location = recode(Location,
                                                                       'CHL' = 'Chiwawa',
                                                                       'CHM' = 'Chumstick',
                                                                       'CHW' = 'Chiwaukum',
                                                                       'ICL' = 'Icicle',
                                                                       'LWN' = 'Little Wenatchee',
                                                                       'MCL' = 'Mission',
                                                                       'NAL' = 'Nason',
                                                                       'PES' = 'Peshastin',
                                                                       'WTL' = 'White River',
                                                                       "GLC" = "Gold",
                                                                       "LBC" = "Libby",
                                                                       "MSH" = "Methow Fish Hatchery",
                                                                       "MRW" = "Upper Methow",
                                                                       "TWR" = "Twisp",
                                                                       "CRW" = "Chewuch",
                                                                       "SCP" = "Spring Creek",
                                                                       "BVC" = "Beaver")) |>
                                              mutate(River = Location) |>
                                              mutate(Type = 'DABOM')) |>
                                  pivot_wider(names_from = Origin,
                                              values_from = c(Spawners,
                                                              Spawners_SE),
                                              names_glue = "{Origin}_{.value}") |>
                                  rlang::set_names(function(x) str_remove(x, "_Spawners")) |>
                                  mutate(Location = factor(Location, levels = c('Below_TUM',
                                                                                'TUM_bb',
                                                                                'Icicle',
                                                                                'Peshastin',
                                                                                'Mission',
                                                                                'Chumstick',
                                                                                'Chiwaukum',
                                                                                'Chiwawa',
                                                                                'Nason',
                                                                                'Little Wenatchee',
                                                                                'White River',
                                                                                "Wenatchee (mainstem)",
                                                                                "Lower Methow",
                                                                                "Gold",
                                                                                "Libby",
                                                                                "Methow Fish Hatchery",
                                                                                "Upper Methow",
                                                                                "Twisp",
                                                                                "Chewuch",
                                                                                "Spring Creek",
                                                                                "Beaver"))) |>
                                  arrange(Location, Type) |>
                                  mutate(h_cv = Hatchery_SE / Hatchery,
                                         n_cv = Natural_SE / Natural,
                                         total_spawners = Natural + Hatchery) |>
                                  select(stream = River,
                                         # Reach,
                                         # reach_type = Type,
                                         total_spawners,
                                         nos_est = Natural,
                                         nos_se = Natural_SE,
                                         nos_cv = n_cv,
                                         hos_est = Hatchery,
                                         hos_se = Hatchery_SE,
                                         hos_cv = h_cv)

                                spwn_yr %<>%
                                  bind_rows(spwn_yr %>%
                                              summarize(stream = "Total",
                                                        across(c(total_spawners,
                                                                 nos_est,
                                                                 hos_est),
                                                               sum),
                                                        across(ends_with("_se"),
                                                               ~ sqrt(sum(.^2, na.rm = T))),
                                                        nos_cv = nos_se / nos_est,
                                                        hos_cv = hos_se / hos_est,
                                                        .groups = "drop")) %>%
                                  rowwise() %>%
                                  mutate(phos = hos_est / (hos_est + nos_est),
                                         phos_se = deltamethod(~ x1 / (x1 + x2),
                                                               mean = c(hos_est,
                                                                        nos_est),
                                                               cov = diag(c(hos_se,
                                                                            nos_se)^2))) %>%
                                  ungroup() %>%
                                  mutate(across(stream,
                                                as.factor),
                                         across(stream,
                                                fct_relevel,
                                                c("Wenatchee (mainstem)",
                                                  "Methow (mainstem)",
                                                  "Total"),
                                                after = Inf)) %>%
                                  arrange(stream)

                                list(results_lst = NA,
                                     fpr = fpr_df,
                                     redd_results = NA,
                                     spwn_rch = NA,
                                     redd_spwn_yr = NA,
                                     spwn_yr = spwn_yr) %>%
                                  return()
                              }

                            },
                            otherwise = NULL))) |>
  select(-file_nm,
         -file_exists)

redd_spwn <- spwn_est %>%
  mutate(redd_spwn_yr = map(results_list,
                            "redd_spwn_yr")) %>%
  select(-results_list) %>%
  unnest(redd_spwn_yr) %>%
  filter(!is.na(observed_redds)) %>%
  makeTableNms()

spwn_df <- spwn_est %>%
  mutate(spwn_yr = map(results_list,
                       "spwn_yr")) %>%
  select(-results_list) %>%
  unnest(spwn_yr) %>%
  makeTableNms()

#-----------------------------------------------------------------
# generate tables based on DABOM output
data("age_table")

dabom_est <- crossing(spawn_year = c(2011:2022)) %>%
  mutate(run_year = spawn_year -1) %>%
  select(run_year, spawn_year) %>%
  mutate(dam_cnt_name = if_else(spawn_year %in% c(2011:2015, 2018),
                                "PriestRapids",
                                "RockIsland")) %>%
  mutate(all_results = map2(spawn_year,
                            dam_cnt_name,
                            .f = function(yr, dam_cnt_name) {

                              cat(paste("\n\t Estimates from", yr, "\n"))

                              # load data
                              load(here("analysis/data/derived_data/estimates",
                                        dam_cnt_name,
                                        paste0("UC_Sthd_DABOM_", yr, ".rda")))

                              # age proportions that year
                              age_prop_yr <- tag_summ |>
                                filter(!group %in% c("Start", "Other")) |>
                                mutate(across(group,
                                              recode,
                                              "BelowPriest" = "Below Priest",
                                              "WellsPool" = "Wells Pool"),
                                       across(origin,
                                              recode,
                                              "H" = "Hatchery",
                                              "W" = "Natural")) |>
                                mutate(run_year = year - 1) |>
                                select(tag_code,
                                       population = group,
                                       origin,
                                       age) |>
                                left_join(age_table |>
                                            select(age,
                                                   fresh_water,
                                                   salt_water,
                                                   total_age),
                                          by = "age") |>
                                pivot_longer(cols = fresh_water:total_age,
                                             names_to = "age_type",
                                             values_to = "age_label") |>
                                group_by(population,
                                         origin,
                                         age_type) |>
                                count(age_label) |>
                                filter(!is.na(age_label)) |>
                                mutate(sample_size = sum(n),
                                       prop = n / sample_size) |>
                                ungroup()

                              # population escapement that year
                              pop_escp_yr <- escape_post %>%
                                filter(param %in% c('LWE',
                                                    'ENL',
                                                    'LMR',
                                                    'OKL', 'FST',
                                                    'ICH', 'JD1', 'JDA', 'PRH', 'PRO', 'PRV', 'RSH', 'TMF')) %>%
                                mutate(param = recode(param,
                                                      'LWE' = 'Wenatchee',
                                                      'ENL' = 'Entiat',
                                                      'LMR' = 'Methow',
                                                      'OKL' = 'Okanogan',
                                                      'FST' = 'Okanogan',
                                                      'ICH' = 'Below Priest',
                                                      'JD1' = 'Below Priest',
                                                      'JDA' = 'Below Priest',
                                                      'PRH' = 'Below Priest',
                                                      'PRO' = 'Below Priest',
                                                      'PRV' = 'Below Priest',
                                                      'RSH' = 'Below Priest',
                                                      'TMF' = 'Below Priest')) %>%
                                mutate(param = factor(param,
                                                      levels = c("Wenatchee",
                                                                 "Entiat",
                                                                 "Methow",
                                                                 "Okanogan",
                                                                 "Below Priest"))) %>%
                                group_by(chain, iter, origin, param) %>%
                                summarize(across(escp,
                                                 sum),
                                          .groups = "drop") %>%
                                group_by(origin,
                                         population = param) %>%
                                summarise(mean = mean(escp),
                                          median = median(escp),
                                          # mode = estMode(escp),
                                          sd = sd(escp),
                                          .groups = 'drop') %>%
                                mutate(across(c(mean, median, sd),
                                              ~ if_else(. < 0, 0, .))) %>%
                                select(-mean) %>%
                                rename(est = median,
                                       se = sd) %>%
                                mutate(across(origin,
                                              recode,
                                              "H" = "hos",
                                              "W" = "nos"),
                                       cv = se / est) %>%
                                pivot_wider(names_from = origin,
                                            values_from = c(est, se, cv),
                                            names_glue = "{origin}_{.value}") %>%
                                mutate(total_escapement = nos_est + hos_est,
                                       phos = hos_est / (hos_est + nos_est)) %>%
                                select(population,
                                       total_escapement,
                                       starts_with("nos"),
                                       starts_with("hos"),
                                       phos)

                              # all escapement locations that year
                              escp_yr = escape_summ %>%
                                mutate(across(origin,
                                              recode,
                                              "H" = "Hatchery",
                                              "W" = "Natural")) %>%
                                select(-species, -spawn_year,
                                       -skew, -kurtosis) %>%
                                mutate(across(mean:mode,
                                              janitor::round_half_up),
                                       across(mean:mode,
                                              as.integer)) %>%
                                rename(estimate = median,
                                       se = sd) %>%
                                select(-mean, -mode)

                              # detection probabilities that year
                              detect_yr <- detect_summ %>%
                                select(node,
                                       n_tags,
                                       estimate = median,
                                       se = sd,
                                       ends_with("CI"))

                              # sex proportions that year
                              sexed_tags <- tag_summ %>%
                                filter(group %in% c("Wenatchee",
                                                    "Entiat",
                                                    "Methow",
                                                    "Okanogan",
                                                    "BelowPriest")) %>%
                                mutate(across(group,
                                              fct_recode,
                                              "Below Priest" = "BelowPriest"),
                                       across(group,
                                              fct_drop)) %>%
                                group_by(population = group,
                                         origin,
                                         sex) %>%
                                summarise(n_tags = n_distinct(tag_code[!is.na(sex)]),
                                          .groups = "drop") %>%
                                filter(!is.na(sex)) %>%
                                full_join(expand(.,
                                                 population,
                                                 origin,
                                                 sex),
                                          by = c("population", "origin", "sex")) %>%
                                mutate(across(n_tags,
                                              replace_na,
                                              0))

                              sex_yr <- sexed_tags %>%
                                bind_rows(sexed_tags %>%
                                            group_by(population, sex) %>%
                                            summarize(across(n_tags,
                                                             sum),
                                                      .groups = "drop") %>%
                                            mutate(origin = "All")) %>%
                                mutate(origin = factor(origin,
                                                       levels = c("W", "H", "All")),
                                       across(origin,
                                              fct_recode,
                                              "Natural" = "W",
                                              "Hatchery" = "H")) %>%
                                mutate(group = paste(population, origin, sep = "_")) %>%
                                group_by(group) %>%
                                mutate(total_sexed = sum(n_tags)) %>%
                                mutate(prop = n_tags / total_sexed,
                                       prop_se = sqrt((prop * (1 - prop)) / total_sexed)) %>%
                                ungroup() %>%
                                arrange(population, origin, sex) %>%
                                left_join(pop_escp_yr %>%
                                            rename(total_est = total_escapement) %>%
                                            rowwise() %>%
                                            mutate(total_se = sqrt(sum(c(nos_se,
                                                                         hos_se)^2))) %>%
                                            ungroup() %>%
                                            select(population,
                                                   ends_with("est"),
                                                   ends_with("se")) %>%
                                            pivot_longer(c(ends_with("est"),
                                                           ends_with("se"))) %>%
                                            mutate(type = if_else(str_detect(name, "_se$"),
                                                                  "se",
                                                                  "est"),
                                                   origin = str_split(name, "_", simplify = T)[,1],
                                                   across(origin,
                                                          recode,
                                                          "total" = "All",
                                                          "nos" = "Natural",
                                                          "hos" = "Hatchery")) %>%
                                            select(-name) %>%
                                            pivot_wider(names_from = type,
                                                        values_from = value),
                                          by = c("population", "origin")) %>%
                                rowwise() %>%
                                mutate(across(est,
                                              ~ . * prop),
                                       se = deltamethod(~ x1 * x2,
                                                        mean = c(est, prop),
                                                        cov = diag(c(se, prop_se)^2))) %>%
                                ungroup() %>%
                                select(population,
                                       origin,
                                       sex,
                                       n_tags,
                                       total_sexed,
                                       starts_with("prop"),
                                       escp = est,
                                       escp_se = se)


                              # convert lengths from mm to cm if needed
                              if(max(tag_summ$fork_length, na.rm = T) > 100) {
                                tag_summ %<>%
                                  mutate(across(fork_length,
                                                measurements::conv_unit,
                                                from = "mm",
                                                to = "cm"))
                              }

                              size_yr <- tag_summ |>
                                filter(!group %in% c("Start", "Other")) |>
                                mutate(across(group,
                                              fct_drop)) |>
                                mutate(across(group,
                                              recode,
                                              "BelowPriest" = "Below Priest",
                                              "WellsPool" = "Wells Pool"),
                                       across(origin,
                                              recode,
                                              "H" = "Hatchery",
                                              "W" = "Natural")) |>
                                select(tag_code,
                                       population = group,
                                       origin,
                                       age,
                                       fork_length) |>
                                left_join(age_table |>
                                            select(age,
                                                   salt_water),
                                          by = "age") |>
                                filter(!is.na(fork_length),
                                       !is.na(salt_water)) |>
                                group_by(population,
                                         origin,
                                         salt_water) |>
                                summarize(across(fork_length,
                                                 list(mean = mean,
                                                      sd = sd,
                                                      n = length),
                                                 .names = "{.fn}"),
                                          .groups = "drop") %>%
                                rename(salt_water_age = salt_water)

                              # proportions / escapement of mark groups, by population
                              mark_tag_summ = tag_summ |>
                                filter(group %in% c("Wenatchee",
                                                    "Entiat",
                                                    "Methow",
                                                    "Okanogan",
                                                    "BelowPriest")) |>
                                mutate(across(group,
                                              fct_drop)) |>
                                mutate(cwt = if_else(is.na(cwt), F,
                                                     if_else(cwt %in% c("SN", "BD"),
                                                             T, NA)),
                                       ad_clip = if_else(is.na(ad_clip) | origin == "W",
                                                         F,
                                                         if_else(ad_clip == "AD", T, NA))) |>
                                mutate(ad_clip_chr = recode(as.character(ad_clip),
                                                            "TRUE" = "AD",
                                                            "FALSE" = "AI"),
                                       cwt_chr = recode(as.character(cwt),
                                                        "TRUE" = "CWT",
                                                        "FALSE" = "noCWT")) |>
                                tidyr::unite("mark_grp", ad_clip_chr, cwt_chr, remove = T) |>
                                mutate(mark_grp = if_else(origin == "W",
                                                          "Wild",
                                                          mark_grp))

                              # proportions of each type of fish (H/W, Ad-clip, CWT combinations) in each population
                              mark_grp_prop = mark_tag_summ |>
                                group_by(population = group,
                                         origin,
                                         ad_clip,
                                         cwt,
                                         mark_grp) |>
                                summarize(n_tags = n_distinct(tag_code),
                                          .groups = "drop") |>
                                full_join(expand(mark_tag_summ,
                                                 population = group,
                                                 nesting(origin,
                                                         ad_clip,
                                                         cwt,
                                                         mark_grp)),
                                          by = c("population", "origin", "ad_clip", "cwt", "mark_grp")) |>
                                arrange(population,
                                        mark_grp) |>
                                mutate(across(n_tags,
                                              replace_na,
                                              0)) |>
                                group_by(population) |>
                                mutate(tot_tags = sum(n_tags),
                                       prop = n_tags / tot_tags,
                                       prop_se = sqrt((prop * (1 - prop))/tot_tags)) |>
                                ungroup() |>
                                arrange(population,
                                        origin,
                                        ad_clip,
                                        cwt)

                              # generate posterior samples of mark proportions
                              set.seed(6)
                              n_iter = max(escape_post$iter)
                              prop_samps = mark_grp_prop |>
                                filter(tot_tags > 0) |>
                                group_by(origin,
                                         population,
                                         tot_tags) |>
                                nest() |>
                                mutate(samp = map(data,
                                                  .f = function(x) {
                                                    rmultinom(n_iter, sum(x$n_tags), x$prop) |>
                                                      set_colnames(1:n_iter) |>
                                                      set_rownames(x$mark_grp) |>
                                                      as_tibble(rownames = 'mark_grp') |>
                                                      pivot_longer(-1,
                                                                   names_to = "iter",
                                                                   values_to = "n_tags") |>
                                                      mutate(across(iter,
                                                                    as.integer)) |>
                                                      group_by(iter) |>
                                                      mutate(prop = n_tags / sum(n_tags)) |>
                                                      arrange(iter, mark_grp) |>
                                                      ungroup()
                                                  })) |>
                                ungroup() |>
                                select(-data) |>
                                unnest(samp) |>
                                left_join(mark_grp_prop |>
                                            select(-n_tags,
                                                   -starts_with("prop")),
                                          by = c("population", "origin", "tot_tags", "mark_grp")) |>
                                bind_rows(mark_grp_prop |>
                                            filter(tot_tags == 0) |>
                                            crossing(iter = 1:n_iter) |>
                                            select(-prop_se)) |>
                                arrange(population, origin, mark_grp, iter) |>
                                select(iter, any_of(names(mark_grp_prop)))

                              # posterior samples
                              mark_post = escape_post |>
                                filter(param %in% c("LWE", 'ENL', 'LMR', 'OKL')) |>
                                mutate(population = recode(param,
                                                           'LWE' = 'Wenatchee',
                                                           'ENL' = 'Entiat',
                                                           'LMR' = 'Methow',
                                                           'OKL' = 'Okanogan')) |>
                                inner_join(prop_samps,
                                           by = c("iter", "origin",
                                                  "population")) |>
                                mutate(across(c(prop,
                                                n_tags,
                                                tot_tags),
                                              replace_na,
                                              0)) |>
                                mutate(prop = if_else(value == 0,
                                                      0,
                                                      prop)) |>
                                mutate(n_fish = escp * prop)

                              mark_grp_yr = mark_post |>
                                group_by(population,
                                         origin,
                                         ad_clip,
                                         cwt,
                                         mark_grp) |>
                                summarise(across(n_fish,
                                                 list(mean = mean,
                                                      median = median,
                                                      # mode = DABOM::estMode,
                                                      se = sd,
                                                      skew = moments::skewness,
                                                      kurtosis = moments::kurtosis,
                                                      lowerCI = ~ coda::HPDinterval(coda::as.mcmc(.x))[,1],
                                                      upperCI = ~ coda::HPDinterval(coda::as.mcmc(.x))[,2]),
                                                 na.rm = T,
                                                 .names = "{.fn}"),
                                          .groups = "drop") |>
                                # mutate(across(mode,
                                #               ~ if_else(. < 0, 0, .))) |>
                                left_join(mark_grp_prop |>
                                            rename(proportion = prop),
                                          by = c("population", "origin", "ad_clip", "cwt", "mark_grp")) |>
                                select(population,
                                       origin:mark_grp,
                                       n_tags:proportion,
                                       prop_se,
                                       everything()) |>
                                mutate(across(population,
                                              factor,
                                              levels = levels(tag_summ$group)),
                                       across(population,
                                              fct_drop)) |>
                                arrange(population, mark_grp) |>
                                mutate(across(origin,
                                              fct_recode,
                                              "Natural" = "W",
                                              "Hatchery" = "H")) |>
                                select(population:prop_se,
                                       estimate = mean,
                                       se,
                                       lowerCI,
                                       upperCI)



                              res_list <- list(age_prop_yr = age_prop_yr,
                                               pop_escp_yr = pop_escp_yr,
                                               escp_yr = escp_yr,
                                               detect_yr = detect_yr,
                                               sex_yr = sex_yr,
                                               size_yr = size_yr,
                                               mark_grp_yr = mark_grp_yr)

                              rm(age_prop_yr,
                                 pop_escp_yr,
                                 escp_yr,
                                 detect_yr,
                                 sexed_tags,
                                 sex_yr,
                                 mark_grp_yr)

                              return(res_list)
                            }))

pop_escp <- dabom_est %>%
  mutate(res = map(all_results,
                   "pop_escp_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) %>%
  makeTableNms()

all_escp <- dabom_est %>%
  mutate(res = map(all_results,
                   "escp_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) %>%
  makeTableNms()

detect_df <- dabom_est %>%
  mutate(res = map(all_results,
                   "detect_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) %>%
  makeTableNms()

sex_prop <- dabom_est %>%
  mutate(res = map(all_results,
                   "sex_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) %>%
  makeTableNms()

size_df <- dabom_est %>%
  mutate(res = map(all_results,
                   "size_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res)

size_tab <- expand(size_df,
                   nesting(run_year, spawn_year),
                   population, origin, salt_water_age) |>
  left_join(size_df,
            by = c("run_year",
                   "spawn_year",
                   "population",
                   "origin",
                   "salt_water_age")) |>
  mutate(across(salt_water_age,
                ~ paste0("Salt-", .))) |>
  rename(Mean = mean,
         SD = sd,
         N = n) |>
  group_by(spawn_year,
           population,
           origin) |>
  mutate(sample_size = sum(N, na.rm = T)) |>
  ungroup() |>
  pivot_wider(names_from = salt_water_age,
              values_from = c(Mean, SD, N),
              names_glue = "{salt_water_age}_{.value}",
              names_vary = "slowest") |>
  makeTableNms()



age_prop <- dabom_est %>%
  mutate(res = map(all_results,
                   "age_prop_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res)

tot_age_prop <- age_prop |>
  filter(age_type == "total_age") |>
  mutate(across(age_label,
                factor,
                levels = 2:11)) |>
  select(-c(age_type, n)) |>
  pivot_wider(names_from = age_label,
              names_expand = T,
              values_from = prop,
              values_fill = 0) |>
  makeTableNms()

sw_age_prop <- age_prop |>
  filter(age_type == "salt_water") |>
  mutate(across(age_label,
                factor,
                levels = 0:4)) |>
  select(-c(age_type, n)) |>
  pivot_wider(names_from = age_label,
              names_expand = T,
              values_from = prop,
              values_fill = 0) |>
  makeTableNms()

fw_age_prop <- age_prop |>
  filter(age_type == "fresh_water") |>
  mutate(across(age_label,
                ~ paste0(., ".X")),
         across(age_label,
                factor,
                levels = paste0(1:5, ".X"))) |>
  select(-c(age_type, n)) |>
  pivot_wider(names_from = age_label,
              names_expand = T,
              values_from = prop,
              values_fill = 0) |>
  makeTableNms()


mark_grp_df <- dabom_est %>%
  mutate(res = map(all_results,
                   "mark_grp_yr")) %>%
  select(-dam_cnt_name,
         -all_results) %>%
  unnest(res) |>
  makeTableNms()


#-----------------------------------------------------------------
# generate table of escapement estimates at Priest
priest_df = tibble(spawn_year = 2011:2022) |>
  mutate(prd_df = map(spawn_year,
                      .f = function(yr) {
                        # load data
                        load(here("analysis/data/derived_data/estimates",
                                  "PriestRapids",
                                  paste0("UC_Sthd_DABOM_", yr, ".rda")))

                        tag_num <- bio_df |>
                          group_by(origin) |>
                          summarize(n_tags = n_distinct(tag_code),
                                    .groups = "drop") |>
                          mutate(total_tags = sum(n_tags))

                        res <- dam_escp_df |>
                          filter(dam == "PriestRapids") |>
                          mutate(run_year = year - 1) |>
                          select(run_year,
                                 # spawn_year = year,
                                 win_cnt,
                                 origin,
                                 contains("reasc_rate"),
                                 contains("priest_cnt")) |>
                          mutate(total = sum(priest_cnt),
                                 total_se = sqrt(sum(priest_cnt_se^2))) |>
                          left_join(tag_num,
                                    by = "origin")
                        return(res)
                      })) |>
  unnest(prd_df) |>
  relocate(run_year,
           .before = 1) |>
  mutate(across(origin,
                recode,
                "H" = "Hatchery",
                "W" = "Natural")) |>
  makeTableNms() |>
  rename(`PRD Ladder Count` = `Win Cnt`,
         `Adult Reascension Adjustment Rate` = `Reasc Rate`,
         `Reascension Rate SE` = `Reasc Rate SE`,
         `Priest Est` = `Priest Cnt`,
         `Priest SE` = `Priest Cnt SE`)


# now a table for years when Rock Island was used
rock_isl_df = priest_df |>
  select(run_year = `Run Year`,
         spawn_year = `Spawn Year`) |>
  filter(!spawn_year %in% c(2011:2015, 2018)) |>
  distinct() |>
  mutate(ria_df = map(spawn_year,
                      .f = function(yr) {
                        # load data
                        load(here("analysis/data/derived_data/estimates",
                                  "RockIsland",
                                  paste0("UC_Sthd_DABOM_", yr, ".rda")))

                        res <- dam_escp_df |>
                          filter(dam == "RockIsland") |>
                          select(origin,
                                 #contains("reasc_rate"),
                                 contains("priest_cnt"),
                                 dam) |>
                          mutate(total = sum(priest_cnt),
                                 total_se = sqrt(sum(priest_cnt_se^2))) |>
                          relocate(dam,
                                   .after = "total_se")
                        return(res)
                      })) |>
  unnest(ria_df) |>
  mutate(across(origin,
                recode,
                "H" = "Hatchery",
                "W" = "Natural"),
         across(dam,
                recode,
                "RockIsland" = "Rock Island")) |>
  makeTableNms() |>
  rename(#`PRD Ladder Count` = `Win Cnt`,
    # `Adult Reascension Adjustment Rate` = `Reasc Rate`,
    # `Reascension Rate SE` = `Reasc Rate_se`,
    `Modelled Priest Est` = `Priest Cnt`,
    `Modelled Priest SE` = `Priest Cnt SE`,
    `Modelled Total` = Total,
    `Modelled Total SE` = `Total SE`,
    `Dam Method` = Dam)


dam_cnt_tab <- priest_df |>
  left_join(rock_isl_df,
            by = c("Run Year",
                   "Spawn Year",
                   "Origin")) |>
  relocate(c(`N Tags`,
             `Total Tags`),
           .after = "Dam Method") |>
  mutate(`Tag Rate` = if_else(is.na(`Dam Method`),
                              `N Tags` / `Priest Est`,
                              `N Tags` / `Modelled Priest Est`),
         `Total Tag Rate` = if_else(is.na(`Dam Method`),
                                    `Total Tags` / Total,
                                    `Total Tags` / `Modelled Total`))


#-----------------------------------------------------------------
# put together in a list to save / write to Excel
# make some decisions about rounding
save_list <- list(
  "Redd Abundance" = redd_spwn |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(c(`Redd Est`,
                    `NOS Est`,
                    `HOS Est`),
                  round_half_up)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Run Escp by Origin" = pop_escp |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(c(`Total Escapement`,
                    ends_with("Est")),
                  round_half_up)) |>
    mutate(`Total Escapement` = `NOS Est` + `HOS Est`) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Spawning Escp by Origin" = spwn_df |>
    bind_rows(pop_escp |>
                filter(Population == "Entiat") |>
                rename(`Total Spawners` = `Total Escapement`) |>
                mutate(Stream = as.factor("Total")) |>
                rowwise() |>
                mutate(`pHOS SE` = deltamethod(~ x1 / (x1 + x2),
                                               mean = c(`HOS Est`,
                                                        `NOS Est`),
                                               cov = diag(c(`HOS SE`,
                                                            `NOS SE`)^2))) |>
                ungroup() |>
                mutate(across(Population,
                              factor,
                              levels = levels(pop_escp$Population)))) |>
    arrange(`Spawn Year`,
            Population,
            Stream) |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(c(`Total Spawners`,
                    ends_with("Est")),
                  round_half_up)) |>
    mutate(`Total Spawners` = `NOS Est` + `HOS Est`) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Run Escp All Locations" = all_escp |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Detection Probability" = detect_df |>
    mutate(across(c(ends_with("Year"),
                    `N Tags`),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Run Escp Sex" = sex_prop |>
    mutate(across(ends_with("Year"),
                  as.integer),
           across(Escp,
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "SaltAge at Maturity" = sw_age_prop |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "FWAge at Maturity" = fw_age_prop |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "TotalAge at Maturity" = tot_age_prop |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Size at Maturity" = size_tab |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Priest Rapids Escapement" = dam_cnt_tab |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(c("PRD Ladder Count",
                    ends_with("Priest Est"),
                    ends_with("Total")),
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3)),

  "Mark Group Population" = mark_grp_df |>
    mutate(across(ends_with("Year"),
                  as.integer)) |>
    mutate(across(Estimate,
                  ~ as.integer(round_half_up(.)))) |>
    mutate(across(where(is.double),
                  round,
                  digits = 3))
)

#-----------------------------------------------------------------
# actually save entire file
write_xlsx(x = save_list,
           path = paste0("T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/Estimates/",
                         "UC_STHD_Model_Output.xlsx"))

#--------------------------------------------------------------
# read in previous estimates, add latest year to them

# what year are we overwriting or adding?
yr = 2022

library(readxl)
output_path <- paste0("T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/Estimates/",
                      "UC_STHD_Model_Output.xlsx")

tab_nms <- excel_sheets(output_path)

previous_output <- as.list(tab_nms) |>
  rlang::set_names() %>%
  map(.f = function(x) {
    read_excel(output_path,
               sheet = x) %>%
      filter(`Spawn Year` != yr)
  })

lst_2022 <- save_list |>
  map(.f = function(x) {
    x |>
      filter(`Spawn Year` == yr)
  })



identical(names(previous_output), names(lst_2022))

save_list = vector("list",
                   length = length(previous_output))
names(save_list) = names(previous_output)
for(i in 1:length(previous_output)) {
  save_list[[i]] <- previous_output[[i]] |>
    bind_rows(lst_2022[[i]])
}

write_xlsx(x = save_list,
           path = paste0("T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/Estimates/",
                         "UC_STHD_Model_Output.xlsx"))

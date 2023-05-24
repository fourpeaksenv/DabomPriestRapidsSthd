# Author: Kevin See
# Purpose: Examine probability of steelhead spawners being detected moving into tributary spawning areas
# Created: 2/9/23
# Last Modified: 2/9/2023
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(janitor)
library(ggrepel)
library(here)

theme_set(theme_bw())

#-----------------------------------------------------------------
# read in estimates of detection probabilities for each node
detect_df <- read_excel("T:/DFW-Team FP Upper Columbia Escapement - General/UC_Sthd/estimates/UC_STHD_Model_Output.xlsx",
                        "Detection Probability") |>
  clean_names()

# estimate detection probability for each site
site_det <- detect_df |>
  mutate(site_code = str_remove(node, "A0$"),
         across(site_code,
                str_remove,
                "B0$")) |>
  # filter(!(str_detect(node, "A0$") | str_detect(node, "B0$")))
  group_by(run_year,
           spawn_year,
           site_code) |>
  mutate(array_nm = str_sub(node, -2),
         array_nm = if_else(!array_nm %in% c("A0", "B0"),
                            "B0",
                            array_nm)) |>
  select(run_year:spawn_year, site_code,
         array_nm,
         estimate) |>
  pivot_wider(names_from = array_nm,
              values_from = estimate) |>
  mutate(across(c(A0, B0),
                replace_na,
                0)) |>
  mutate(p_hat = 1 - ((1 - A0) * (1 - B0))) |>
  filter(p_hat > 0,
         p_hat < 1) |>
  group_by(site_code) |>
  summarize(across(p_hat,
                   list(mean = mean,
                        median = median)))

# group sites by subbasin
plot_df <- site_det |>
  filter(site_code %in% c("ENL",
                          "ENA",
                          "MAD")) |>
  mutate(across(contains("p_hat"),
                ~ if_else(site_code == "ENL",
                          .,
                          1 - (1 - .[site_code == "ENL"]) * (1 - .)))) |>
  arrange(p_hat_median) |>
  mutate(subbasin = "Entiat") |>
  bind_rows(site_det |>
              filter(site_code %in% c("LWE",
                                      "PES",
                                      "MCL",
                                      "CHM",
                                      "TUM",
                                      "ICL")) |>
              mutate(across(contains("p_hat"),
                            ~ if_else(site_code == "LWE",
                                      .,
                                      1 - (1 - .[site_code == "LWE"]) * (1 - .)))) |>

              arrange(p_hat_median) |>
              mutate(subbasin = "Wenatchee")) |>
  bind_rows(site_det |>
              filter(site_code %in% c("LMR",
                                      "GLC",
                                      "LBC",
                                      # "MRC",
                                      "MSH",
                                      "MRW",
                                      "TWR",
                                      "CRW",
                                      "SCP",
                                      "BVC")) |>
              mutate(across(contains("p_hat"),
                            ~ if_else(site_code == "LMR",
                                      .,
                                      1 - (1 - .[site_code == "LMR"]) * (1 - .)))) |>
              arrange(p_hat_median) |>
              mutate(subbasin = "Methow")) |>
  bind_rows(site_det |>
              filter(site_code %in% c("OKL",
                                      "AEN",
                                      "WAN",
                                      "LLC",
                                      "TNK",
                                      "ANT",
                                      "BPC",
                                      "WHS",
                                      "JOH",
                                      "ZSL",
                                      "OMK",
                                      "SA1")) |>
              mutate(across(contains("p_hat"),
                            ~ if_else(site_code == "OKL",
                                      .,
                                      1 - (1 - .[site_code == "OKL"]) * (1 - .)))) |>
              arrange(p_hat_median) |>
              mutate(subbasin = "Okanogan"))

ggplot(plot_df,
       aes(x = subbasin,
           y = p_hat_mean,
           fill = subbasin)) +
  geom_boxplot() +
  labs(x = "Population",
       y = "Avg. Probability of Detection") +
  scale_fill_brewer(palette = "Set1",
                    name = "Subbasin") +
  theme(legend.position = "none") +
  geom_label_repel(data = plot_df |>
                     filter(site_code %in% c("ENL",
                                             "LWE",
                                             "LMR",
                                             "OKL")),
                   aes(label = site_code),
                   fill = "lightgray")

ggsave(here("analysis/figures",
            "spawner_detection.png"),
       width = 8,
       height = 5)

plot_df |>
  select(subbasin,
         site_code,
         det_prob = p_hat_mean) |>
  write_csv(here("analysis/data/derived_data",
                 "prob_detection.csv"))

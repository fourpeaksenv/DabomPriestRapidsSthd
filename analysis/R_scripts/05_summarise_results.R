# Author: Kevin See
# Purpose: summarise DABOM results
# Created: 4/1/20
# Last Modified: 4/1/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(DABOM)
library(PITcleanr)
library(tidyverse)
# library(jagsUI)
library(STADEM)
library(readxl)
library(WriteXLS)
library(msm)
library(moments)
library(coda)

#-----------------------------------------------------------------
# set species
spp = "Chinook"
# set year
yr = 2019

#-----------------------------------------------------------------
# load JAGS MCMC results
load(paste0("analysis/data/derived_data/model_fits/TUM_DABOM_", spp, '_', yr,'.rda'))

bio_df = read_rds('analysis/data/derived_data/Bio_Data_2008_2019.rds') %>%
  filter(Year == yr,
         TagID %in% unique(proc_list$ProcCapHist$TagID))

# estimate final spawning location
tag_summ = summariseTagData(proc_list$ProcCapHist,
                            trap_data = bio_df)

# summarise detection probabilities
detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   capHist_proc = proc_list$proc_ch) %>%
  filter(!is.na(mean))

# which sites had detection probabilities fixed at 0% or 100%
detect_summ %>%
  filter(sd == 0)

# look at all the other sites
detect_summ %>%
  filter(sd > 0) %>%
  arrange(desc(sd))

# compile all movement probabilities, and multiply them appropriately
trans_df = compileTransProbs_TUM(dabom_mod)

# summarize transition probabilities
trans_summ = trans_df %>%
  group_by(Origin, param) %>%
  summarise(mean = mean(value),
            median = median(value),
            mode = estMode(value),
            sd = sd(value),
            skew = moments::skewness(value),
            kurtosis = moments::kurtosis(value),
            lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2]) %>%
  mutate_at(vars(mean, median, mode, sd, matches('CI$')),
            list(~ if_else(. < 0, 0, .))) %>%
  ungroup()

#-----------------------------------------------------------------
# total escapement past Tumwater
escp = read_excel(paste0('analysis/data/raw_data/WDFW/PIT_Cleaned_', yr, ' SPCH Data.xlsx'),
                  3,
                  range = anchored("C6",
                                   dim = c(3, 15)),
                  col_names = c("Fish",
                                paste(rep('H', 7), c(rep('M', 4), rep('F', 3)), c(2:5, 3:5), sep = '_'),
                                paste(rep('W', 6), c(rep('M', 3), rep('F', 3)), c(3:5, 3:5), sep = '_'),
                                "Totals"))

# we know exact escapement for hatchery and wild fish. Tags are not proportional to hatchery/wild, so don't use those.
org_escape = escp %>%
  select(-Totals) %>%
  gather(type, nFish, -Fish) %>%
  filter(!is.na(Fish)) %>%
  mutate(Origin = if_else(grepl('H', type),
                          'Hatchery', 'Natural')) %>%
  group_by(Origin) %>%
  summarise(tot_escp = sum(nFish),
            tot_escp_se = 0) %>%
  mutate(Species = 'Chinook',
         SpawnYear = yr) %>%
  select(Species, SpawnYear, everything())

# translate movement estimates to escapement
escape_summ = trans_df %>%
  left_join(org_escape %>%
              select(Origin, tot_escp)) %>%
  mutate(escp = value * tot_escp) %>%
  group_by(Origin, location = param) %>%
  summarise(mean = mean(escp),
            median = median(escp),
            mode = estMode(escp),
            sd = sd(escp),
            skew = moments::skewness(escp),
            kurtosis = moments::kurtosis(escp),
            lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
            upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2]) %>%
  mutate_at(vars(mean, median, mode, sd, matches('CI$')),
            list(~ if_else(. < 0, 0, .))) %>%
  ungroup() %>%
  mutate(Species = spp,
         SpawnYear = yr) %>%
  select(Species, SpawnYear, everything())


bioSumm = tag_summ %>%
  group_by(Stream = Group,
           Origin,
           Sex,
           Age) %>%
  summarise(nTags = n_distinct(TagID),
            meanPOH = mean(POH, na.rm = T)) %>%
  ungroup() %>%
  full_join(expand.grid(list(Stream = levels(tag_summ$Group),
                             Origin = unique(tag_summ$Origin),
                             Sex = unique(tag_summ$Sex),
                             Age = unique(tag_summ$Age)))) %>%
  mutate_at(vars(nTags),
            funs(if_else(is.na(.),
                         as.integer(0), .))) %>%
  group_by(Stream, Origin) %>%
  mutate(totTags = sum(nTags)) %>%
  ungroup() %>%
  mutate(prop = nTags / totTags,
         propSE = sqrt((prop * (1 - prop)) / (totTags))) %>%
  mutate(Stream = as.character(Stream)) %>%
  select(Stream, Origin, totTags, Sex, Age, nTags, prop, propSE, meanPOH) %>%
  arrange(Stream, Origin, Sex, Age)

fullSumm = escape_summ %>%
  filter(location %in% c('past_CHL',
                     'past_CHW',
                     'past_NAL',
                     'past_LWN',
                     'past_WTL',
                     'TUM_bb')) %>%
  mutate(Stream = recode(location,
                         'past_CHL' = 'Chiwawa',
                         'past_CHW' = 'Chiwaukum',
                         'past_NAL' = 'Nason',
                         'past_LWN' = 'LittleWenatchee',
                         'past_WTL' = 'WhiteRiver')) %>%
  mutate(Origin = recode(Origin,
                         'Natural' = 'W',
                         'Hatchery' = 'H')) %>%
  select(Species, SpawnYear, Stream, Origin, Escape = mean, EscSE = sd) %>%
  left_join(bioSumm) %>%
  rowwise() %>%
  mutate(Est = Escape * prop,
         EstSE = deltamethod(~ x1 * x2,
                             mean = c(Escape, prop),
                             cov = diag(c(EscSE, propSE)^2))) %>%
  ungroup() %>%
  mutate(Est = round(Est),
         Est = as.integer(Est)) %>%
  select(Species:Origin, Sex, Age, nTags, meanPOH, prop, Est, EstSE)

# sex proportions by origin
sexOrigSumm = escape_summ %>%
  filter(location %in% c('past_CHL',
                     'past_CHW',
                     'past_NAL',
                     'past_LWN',
                     'past_WTL',
                     'TUM_bb')) %>%
  mutate(Stream = recode(location,
                         'past_CHL' = 'Chiwawa',
                         'past_CHW' = 'Chiwaukum',
                         'past_NAL' = 'Nason',
                         'past_LWN' = 'LittleWenatchee',
                         'past_WTL' = 'WhiteRiver')) %>%
  mutate(Origin = recode(Origin,
                         'Natural' = 'W',
                         'Hatchery' = 'H')) %>%
  select(Species, SpawnYear, Stream, Origin, Escape = mean, EscSE = sd) %>%
  left_join(tag_summ %>%
              filter(!is.na(Group)) %>%
              group_by(Group, Origin, Sex) %>%
              summarise(nTags = n_distinct(TagID[!is.na(Sex)])) %>%
              filter(!is.na(Sex)) %>%
              ungroup() %>%
              spread(Sex, nTags,
                     fill = as.integer(0)) %>%
              mutate(total_sexed = F + M,
                     propF = F / (F + M),
                     propM = 1 - propF,
                     propSexSE = sqrt((propF * (1 - propF)) / (M + F))) %>%
              select(Group, Origin, total_sexed:propSexSE) %>%
              rename(Stream = Group) %>%
              gather(Sex, prop, propF:propM) %>%
              mutate(Sex = recode(Sex,
                                  'propF' = 'F',
                                  'propM' = 'M')) %>%
              select(Stream, Origin, Sex, total_sexed, prop, propSexSE) %>%
              arrange(Stream, Origin, Sex)) %>%
  rowwise() %>%
  mutate(Est = Escape * prop,
         Est_SE = deltamethod(~ x1 * x2,
                              mean = c(Escape, prop),
                              cov = diag(c(EscSE, propSexSE)^2))) %>%
  ungroup() %>%
  mutate(total_sexed = if_else(is.na(total_sexed),
                               as.integer(0), total_sexed)) %>%
  select(Species:Origin, Sex, nSexed = total_sexed, matches('^prop'), matches('Est')) %>%
  arrange(Species, SpawnYear, Stream, Origin, Sex)

sexOrigSumm %<>%
  bind_rows(sexOrigSumm %>%
              group_by(Species, SpawnYear, Stream, Sex) %>%
              summarise(nSexed = sum(nSexed),
                        Est = sum(round(Est)),
                        Est_SE = sqrt(sum(Est_SE^2))) %>%
              ungroup() %>%
              mutate(Origin = 'All')) %>%
  mutate(Origin = factor(Origin,
                         levels = c('W', 'H', 'All'))) %>%
  arrange(Species, SpawnYear, Stream, Origin, Sex) %>%
  mutate(Est = as.integer(round(Est)),
         Est_SE = round(Est_SE, 1))

# generate population level estimates
pop_summ = escape_summ %>%
  filter(location %in% c('past_CHL',
                         'past_CHW',
                         'past_NAL',
                         'past_LWN',
                         'past_WTL',
                         'TUM_bb')) %>%
  mutate(Stream = recode(location,
                         'past_CHL' = 'Chiwawa',
                         'past_CHW' = 'Chiwaukum',
                         'past_NAL' = 'Nason',
                         'past_LWN' = 'LittleWenatchee',
                         'past_WTL' = 'WhiteRiver')) %>%
  mutate(Origin = recode(Origin,
                         'Natural' = 'W',
                         'Hatchery' = 'H')) %>%
  select(Species, SpawnYear, Stream, everything(), -location)

# biological summaries, based on tags detected within each stream / population
bio_list = list('Origin' = tag_summ %>%
                 filter(!is.na(Group)) %>%
                 group_by(Group, Origin) %>%
                 summarise(nTags = n_distinct(TagID)) %>%
                 spread(Origin, nTags,
                        fill = as.integer(0)) %>%
                 ungroup() %>%
                 mutate(propW = W / (W + H),
                        propH = 1 - propW,
                        propOrgSE = sqrt((propW * (1 - propW)) / (W + H))),
               'AllSex' = tag_summ %>%
                 filter(!is.na(Group)) %>%
                 group_by(Group, Origin, Sex) %>%
                 summarise(nTags = n_distinct(TagID[!is.na(Sex)])) %>%
                 filter(!is.na(Sex)) %>%
                 ungroup() %>%
                 spread(Sex, nTags,
                        fill = as.integer(0)) %>%
                 mutate(total_sexed = F + M,
                        propF = F / (F + M),
                        propM = 1 - propF,
                        propSexSE = sqrt((propF * (1 - propF)) / (M + F))) %>%
                 select(Group, Origin, total_sexed:propSexSE),
               'AllAge' = tag_summ %>%
                 filter(!is.na(Group)) %>%
                 group_by(Group, Origin, Age) %>%
                 summarise(nTags = n_distinct(TagID)) %>%
                 ungroup() %>%
                 group_by(Group, Origin) %>%
                 mutate(total_aged = sum(nTags)) %>%
                 spread(Age, nTags,
                        fill = 0) %>%
                 select(Group, Origin, total_aged,
                        # not_aged = `<NA>`,
                        everything()) %>%
                 mutate_at(vars(total_aged:`5`),
                           funs(as.integer)))

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
                     mutate_at(vars(-Node, -n_tags),
                               list(round),
                               digits = 3) %>%
                     rename(estimate = mean,
                            se = sd) %>%
                     select(-median, -mode),
                   'Sex by Origin by Stream' = sexOrigSumm,
                   'Biological Summary' = fullSumm),
              bio_list)

WriteXLS(x = save_list,
         ExcelFileName = paste0('outgoing/estimates/TUM_Chnk_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx'),
         AdjWidth = T,
         AutoFilter = F,
         BoldHeaderRow = T,
         FreezeRow = 1)

# Author: Kevin See
# Purpose: summarise DABOM results
# Created: 4/1/20
# Last Modified: 10/29/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(DABOM)
library(PITcleanr)
library(tidyverse)
# library(jagsUI)
library(STADEM)
library(readxl)
# library(WriteXLS)
library(openxlsx)
library(msm)
library(moments)
library(coda)

#-----------------------------------------------------------------
# set species
spp = "Steelhead"
# set year
yr = 2019

for(yr in 2011:2019) {
  #-----------------------------------------------------------------
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

  # summarise detection probabilities
  detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                     capHist_proc = proc_list$proc_ch %>%
                                       filter(UserProcStatus)) %>%
    filter(!is.na(mean))


  # which sites had detection probabilities fixed at 0% or 100%
  detect_summ %>%
    filter(sd == 0) %>%
    arrange(Node, mean, n_tags)

  # look at all the other sites
  detect_summ %>%
    filter(sd > 0) %>%
    arrange(desc(sd))

  # compile all movement probabilities, and multiply them appropriately
  trans_df = compileTransProbs_PRA(dabom_mod)

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
  # total escapement past Priest
  #-----------------------------------------------------------------
  # window count
  tot_win_cnt = getWindowCounts(dam = 'PRD',
                                spp = spp,
                                start_date = paste0(yr-1, '0601'),
                                end_date = paste0(yr, '0531')) %>%
    summarise_at(vars(win_cnt),
                 list(sum)) %>%
    pull(win_cnt)

  # reascension data
  reasc_data = queryPITtagData(damPIT = 'PRA',
                               spp = spp,
                               start_date = paste0(yr-1, '0601'),
                               end_date = paste0(yr, '0531'))

  # adjust for re-ascension and origin
  tot_escape = reasc_data %>%
    mutate(SpawnYear = yr,
           TagIDAscentCount = ifelse(is.na(TagIDAscentCount),
                                     0, TagIDAscentCount),
           ReAscent = ifelse(TagIDAscentCount > 1, T, F)) %>%
    group_by(Species, SpawnYear, Date) %>%
    summarise(tot_tags = n_distinct(TagID),
              reascent_tags = n_distinct(TagID[ReAscent])) %>%
    ungroup() %>%
    group_by(Species, SpawnYear) %>%
    summarise_at(vars(matches('tags')),
                 list(sum),
                 na.rm = T) %>%
    ungroup() %>%
    mutate(reascRate = reascent_tags / tot_tags,
           reascRateSE = sqrt(reascRate * (1 - reascRate) / tot_tags),
           totWinCnt = tot_win_cnt,
           adjWinCnt = tot_win_cnt * (1 - reascRate),
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
                       propOrgSE = sqrt((propW * (1 - propW)) / (W + H))))

  org_escape = tot_escape %>%
    mutate(Hescp = propH * adjWinCnt,
           HescpSE = deltamethod(~ x1 * x2,
                                 mean = c(propH, adjWinCnt),
                                 cov = diag(c(propOrgSE, adjWinCntSE)^2))) %>%
    mutate(Wescp = propW * adjWinCnt,
           WescpSE = deltamethod(~ x1 * x2,
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



  #-----------------------------------------------------------------
  # translate movement estimates to escapement
  #-----------------------------------------------------------------
  # bootstrap posteriors
  n_samps = trans_df %>%
    group_by(Origin, param) %>%
    summarise(nIters = n()) %>%
    ungroup() %>%
    pull(nIters) %>%
    unique()

  set.seed(5)
  escape_summ = org_escape %>%
    split(list(.$Origin)) %>%
    map_df(.id = 'Origin',
           .f = function(x) {
             tibble(totEsc = rnorm(n_samps, x$tot_escp, x$tot_escp_se)) %>%
               mutate(iter = 1:n_samps)
           }) %>%
    left_join(trans_df %>%
                group_by(Origin, param) %>%
                mutate(iter = 1:n()) %>%
                select(-chain) %>%
                rename(prob = value) %>%
                slice(sample.int(max(iter), n_samps)) %>%
                ungroup()) %>%
    mutate(value = totEsc * prob) %>%
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
              funs(ifelse(. < 0, 0, .))) %>%
    ungroup()


  # generate population level estimates
  pop_summ = escape_summ %>%
    filter(param %in% c("past_LWE", 'past_ENL', 'past_LMR', 'past_OKL', 'dwnStrm', 'WEA_bb')) %>%
    mutate(Population = recode(param,
                               'past_LWE' = 'Wenatchee',
                               'past_ENL' = 'Entiat',
                               'past_LMR' = 'Methow',
                               'past_OKL' = 'Okanogan',
                               'dwnStrm' = 'BelowPriest',
                               'WEA_bb' = 'WellsPool')) %>%
    arrange(Origin, Population) %>%
    select(Origin, Population, everything(), -param)

  #-----------------------------------------------------------------
  # combine biological data with escapement estimates
  #-----------------------------------------------------------------

  bioSumm = tag_summ %>%
    group_by(Stream = Group,
             Origin,
             Sex,
             Age) %>%
    summarise(nTags = n_distinct(TagID),
              meanFL = mean(ForkLength, na.rm = T)) %>%
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
    select(Stream, Origin, totTags, Sex, Age, nTags, prop, propSE, meanFL) %>%
    arrange(desc(Stream), Origin, Sex, Age)

  fullSumm = pop_summ %>%
    mutate(Origin = recode(Origin,
                           'Natural' = 'W',
                           'Hatchery' = 'H')) %>%
    mutate(Species = spp,
           SpawnYear = yr) %>%
    select(Species, SpawnYear, Population, Origin,
           Escape = mean,
           EscSE = sd) %>%
    left_join(bioSumm) %>%
    rowwise() %>%
    mutate(Est = Escape * prop,
           EstSE = deltamethod(~ x1 * x2,
                               mean = c(Escape, prop),
                               cov = diag(c(EscSE, propSE)^2))) %>%
    ungroup() %>%
    mutate(Est = round(Est),
           Est = as.integer(Est)) %>%
    select(Species, SpawnYear, Population, Origin, Sex, Age, nTags, meanFL, prop, Est, EstSE)

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
                    mutate_at(vars(-(Group:Origin)),
                              list(as.integer)))

  #-----------------------------------------------------------------
  # write results to an Excel file
  save_list = c(list('Population Escapement' = pop_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate_at(vars(mean:mode),
                                 list(round),
                                 digits = 0) %>%
                       mutate_at(vars(sd:upperCI),
                                 list(round),
                                 digits = 1) %>%
                       rename(estimate = mean,
                              se = sd) %>%
                       select(-median, -mode),
                     'All Escapement' = escape_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate_at(vars(mean:mode),
                                 list(round),
                                 digits = 0) %>%
                       mutate_at(vars(sd:upperCI),
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
                     'Biological Summary' = fullSumm),
                bio_list)

  # WriteXLS(x = save_list,
  #          ExcelFileName = paste0('outgoing/estimates/PRA_', spp, '_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx'),
  #          AdjWidth = T,
  #          AutoFilter = F,
  #          BoldHeaderRow = T,
  #          FreezeRow = 1)

  # using a different package to write to Excel
  openxlsx::write.xlsx(x = save_list,
                       file = paste0('outgoing/estimates/PRA_', spp, '_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx'),
                       firstRow = T,
                       headerStyle = openxlsx::createStyle(textDecoration = "BOLD"))
}

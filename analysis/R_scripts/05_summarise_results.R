# Author: Kevin See
# Purpose: summarise DABOM results
# Created: 4/1/20
# Last Modified: 1/12/21
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(DABOM)
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

# for years prior to 2020, use older version of DABOM
if(yr < 2020) {
  remotes::install_github("BiomarkABS/DABOM@v1.0.0")
  library(DABOM)
} else {
  remotes::install_github("BiomarkABS/DABOM")
  library(DABOM)
}


# for(yr in 2011:2019) {
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
                                slice(1),
                              saveCSV = F,
                              file_name = paste0('outgoing/other/TagSummary_', yr, '.csv')) %>%
    mutate(Group = fct_explicit_na(Group),
           Group = fct_relevel(Group,
                               "BelowPriest",
                               after = 4),
           Age = factor(Age,
                        levels = sort(unique(Age))))


  # any duplicate tags?
  sum(duplicated(tag_summ$TagID))
  tag_summ %>%
    filter(TagID %in% TagID[duplicated(TagID)]) %>%
    arrange(TagID, TrapDate)

  # summarize detection probabilities
  detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                     capHist_proc = proc_list$proc_ch %>%
                                       filter(UserProcStatus)) %>%
    filter(!is.na(mean)) %>%
    rename(se = sd)


  # which sites had detection probabilities fixed at 0% or 100%
  detect_summ %>%
    filter(se == 0) %>%
    arrange(Node, mean, n_tags)

  # look at all the other sites
  detect_summ %>%
    filter(se > 0) %>%
    arrange(desc(se))

  # compile all movement probabilities, and multiply them appropriately
  trans_df = compileTransProbs_PRA(dabom_mod)

  # summarize transition probabilities
  trans_summ = trans_df %>%
    group_by(Origin, param) %>%
    summarise(mean = mean(value),
              median = median(value),
              mode = estMode(value),
              se = sd(value),
              skew = moments::skewness(value),
              kurtosis = moments::kurtosis(value),
              lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
              upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2],
              .groups = "drop") %>%
    mutate(across(c(mean, median, mode, se, matches('CI$')),
              ~ if_else(. < 0, 0, .)))

  #-----------------------------------------------------------------
  # total escapement past Priest
  #-----------------------------------------------------------------
  # get escapement past Priest by origin
  start_date = paste0(yr-1, '0601')
  end_date = paste0(yr, '0531')


  # start with PIT-tag based reascension data
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
    # add window counts
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
                            values_fill = list(nTags = as.integer(0)),
                            names_sort = T) %>%
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

  #-----------------------------------------------------------------
  # translate movement estimates to escapement
  #-----------------------------------------------------------------
  # bootstrap posteriors
  n_samps = trans_df %>%
    group_by(Origin, param) %>%
    summarise(nIters = n(),
              .groups = "drop") %>%
    pull(nIters) %>%
    unique()

  if(length(n_samps) > 1) {
    cat("Not every parameter has the same number of posterior samples.\n")
  }

  set.seed(5)
  escp_post = org_escape %>%
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
                 rename(prob = value)) %>%
    mutate(value = prob * escp) %>%
    mutate(Site = if_else(grepl('^past', param),
                          str_remove(param, "past_"),
                          NA_character_))

  escape_summ = escp_post %>%
    group_by(Origin, param) %>%
    summarise(mean = mean(value),
              median = median(value),
              mode = estMode(value),
              se = sd(value),
              skew = moments::skewness(value),
              kurtosis = moments::kurtosis(value),
              lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
              upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2],
              .groups = 'drop') %>%
    mutate(across(c(mean, median, mode, se, matches('CI$')),
                  ~ifelse(. < 0, 0, .)))

  # generate population level estimates
  pop_summ = escape_summ %>%
    filter(param %in% c("past_LWE", 'past_ENL', 'past_LMR', 'past_OKL', 'dwnStrm', 'WEA_bb')) %>%
    mutate(Population = recode(param,
                               'past_LWE' = 'Wenatchee',
                               'past_ENL' = 'Entiat',
                               'past_LMR' = 'Methow',
                               'past_OKL' = 'Okanogan',
                               'dwnStrm' = 'BelowPriest',
                               'WEA_bb' = 'WellsPool'),
           Population = factor(Population,
                               levels = c('Wenatchee',
                                          'Entiat',
                                          'Methow',
                                          'Okanogan',
                                          'BelowPriest',
                                          'WellsPool'))) %>%
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
           EscSE = se) %>%
    left_join(bioSumm %>%
                rename(Population = Stream)) %>%
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
  bio_list = list('All Origin' = tag_summ %>%
                    filter(!is.na(Group)) %>%
                    group_by(Group, Origin) %>%
                    summarise(nTags = n_distinct(TagID)) %>%
                    spread(Origin, nTags,
                           fill = as.integer(0)) %>%
                    ungroup() %>%
                    mutate(propW = W / (W + H),
                           propH = 1 - propW,
                           propOrgSE = sqrt((propW * (1 - propW)) / (W + H))),
                  'All Sex' = tag_summ %>%
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
                           propSexSE = sqrt((propF * (1 - propF)) / (M + F))),
                  'All Age' = tag_summ %>%
                    filter(!is.na(Group)) %>%
                    group_by(Group, Origin, Age) %>%
                    summarise(nTags = n_distinct(TagID),
                              .groups = "drop") %>%
                    group_by(Group, Origin) %>%
                    mutate(total_aged = sum(nTags)) %>%
                    ungroup() %>%
                    pivot_wider(names_from = "Age",
                                values_from = "nTags",
                                values_fill = 0,
                                names_sort = T) %>%
                    select(Group, Origin, total_aged,
                           # not_aged = `<NA>`,
                           everything()) %>%
                    mutate_at(vars(-(Group:Origin)),
                              list(as.integer)))

  #-----------------------------------------------------------------
  # break down hatchery estimates by mark type
  #-----------------------------------------------------------------
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

  # biological summaries of each mark group
  mrk_tag_summ = tag_summ %>%
    filter(Origin == "H") %>%
    mutate(CWT = if_else(is.na(CWT), F,
                         if_else(CWT %in% c("SN", "BD"),
                                 T, NA)),
           AdClip = if_else(is.na(AdClip), F,
                            if_else(AdClip == "AD", T, NA))) %>%
    mutate(AdClip_chr = recode(as.character(AdClip),
                               "TRUE" = "AD",
                               "FALSE" = "AI"),
           CWT_chr = recode(as.character(CWT),
                            "TRUE" = "CWT",
                            "FALSE" = "noCWT")) %>%
    tidyr::unite("Mark_Grp", AdClip_chr, CWT_chr, remove = T)

  mrk_bioSumm = mrk_tag_summ %>%
    group_by(Stream = Group,
             Origin,
             CWT,
             AdClip,
             Mark_Grp,
             Sex,
             Age) %>%
    summarise(nTags = n_distinct(TagID),
              meanFL = mean(ForkLength, na.rm = T),
              .groups = "drop") %>%
    # full_join(expand(mrk_tag_summ,
    #                  Stream = Group,
    #                  Origin,
    #                  CWT,
    #                  AdClip,
    #                  Sex,
    #                  Age)) %>%
    mutate(across(nTags,
                  replace_na,
                  0)) %>%
    group_by(Stream,
             Origin,
             CWT,
             AdClip,
             Mark_Grp) %>%
    mutate(totTags = sum(nTags)) %>%
    ungroup() %>%
    mutate(prop = nTags / totTags,
           propSE = sqrt((prop * (1 - prop)) / (totTags))) %>%
    mutate(Stream = as.character(Stream)) %>%
    select(Stream, Origin, CWT, AdClip, Mark_Grp, Sex, Age, totTags, nTags, prop, propSE, meanFL) %>%
    arrange(desc(Stream), Origin, CWT, AdClip, Sex, Age)

  # biological summaries, based on tags detected within each stream / population
  mrk_bio_list = list('Mark Group Bio Summary' = mrk_bioSumm,
                      'Mark Group Sex' = mrk_tag_summ %>%
                        filter(!is.na(Group)) %>%
                        group_by(Group,
                                 AdClip,
                                 CWT,
                                 Mark_Grp,
                                 Sex) %>%
                        summarise(nTags = n_distinct(TagID[!is.na(Sex)]),
                                  .groups = "drop") %>%
                        filter(!is.na(Sex)) %>%
                        pivot_wider(names_from = "Sex",
                                    values_from = "nTags",
                                    values_fill = 0,
                                    names_sort = T) %>%
                        mutate(total_sexed = F + M,
                               propF = F / (F + M),
                               propM = 1 - propF,
                               propSexSE = sqrt((propF * (1 - propF)) / (M + F))),
                      'Mark Group Age' = mrk_tag_summ %>%
                        filter(!is.na(Group)) %>%
                        group_by(Group,
                                 AdClip,
                                 CWT,
                                 Mark_Grp,
                                 Age) %>%
                        summarise(nTags = n_distinct(TagID),
                                  .groups = "drop") %>%
                        group_by(Group, Mark_Grp) %>%
                        mutate(total_aged = sum(nTags)) %>%
                        ungroup() %>%
                        pivot_wider(names_from = "Age",
                                    values_from = "nTags",
                                    values_fill = 0,
                                    names_sort = T) %>%
                        select(Group:Mark_Grp, total_aged,
                               # not_aged = `<NA>`,
                               everything()) %>%
                        mutate(across(-c(Group:Mark_Grp),
                                      as.integer)))



  # proportions of each type of fish (Ad-clip, CWT combinations) past each site
  mark_prop = mrk_tag_summ %>%
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
    group_by(path_sites, AdClip, CWT, Mark_Grp) %>%
    summarise(n_tags = n_distinct(TagID),
              .groups = "drop") %>%
    right_join(crossing(path_sites = valid_paths$Site,
                        expand(mrk_tag_summ, nesting(AdClip, CWT, Mark_Grp)))) %>%
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

  # generate posterior samples of mark proportions
  set.seed(6)
  prop_samps = mark_prop %>%
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
    mutate(Origin = 'Hatchery') %>%
    select(Origin,
           Site:CWT,
           Mark_Grp,
           n_tags, tot_tags,
           chain, samps) %>%
    unnest(cols = samps)

  # generate escapement estimates by mark group
  # posterior samples
  mark_post = escp_post %>%
    filter(!is.na(Site)) %>%
    filter(Origin == "Hatchery") %>%
    left_join(prop_samps) %>%
    mutate(across(c(prop,
                    n_tags,
                    tot_tags),
                  replace_na,
                  0)) %>%
    mutate(prop = if_else(Origin != "Hatchery",
                          1,
                          if_else(value == 0,
                                  0,
                                  prop))) %>%
    mutate(n_fish = value * prop) %>%
    mutate(across(n_fish,
                  round,
                  digits = 2))

  mark_grps = mark_post %>%
    group_by(Species,
             Origin,
             Site,
             AdClip,
             CWT,
             Mark_Grp,
             n_tags,
             tot_tags) %>%
    summarise(across(n_fish,
                     tibble::lst(mean,
                                 median,
                                 mode = estMode,
                                 se = sd,
                                 skew = moments::skewness,
                                 kurtosis = moments::kurtosis),
                     na.rm = T),
              .groups = "drop") %>%
    mutate(across(n_fish_mode,
                  ~ if_else(. < 0, 0, .))) %>%
    full_join(mark_post %>%
                group_by(Species,
                         Origin,
                         Site,
                         AdClip,
                         CWT,
                         Mark_Grp) %>%
                summarise(lowerCI = coda::HPDinterval(coda::as.mcmc(n_fish))[,1],
                          upperCI = coda::HPDinterval(coda::as.mcmc(n_fish))[,2],
                          .groups = "drop")) %>%
    # select(Species:CWT,
    #        Mark_Grp,
    #        everything()) %>%
    tibble::add_column(proportion = NA_real_,
                       .after = "tot_tags") %>%
    mutate(proportion = n_tags / tot_tags)

  names(mark_grps) = str_remove(names(mark_grps), "^n_fish_")

  pop_mark_grp = mark_grps %>%
    filter(Site %in% c("LWE", 'ENL', 'LMR', 'OKL')) %>%
    mutate(Population = recode(Site,
                               'LWE' = 'Wenatchee',
                               'ENL' = 'Entiat',
                               'LMR' = 'Methow',
                               'OKL' = 'Okanogan'),
           Population = factor(Population,
                               levels = c('Wenatchee',
                                          'Entiat',
                                          'Methow',
                                          'Okanogan'))) %>%
    arrange(Origin, Population) %>%
    select(Origin, Population, everything(), -Site)

  mark_grp_list = c(list('Mark Group Population' = pop_mark_grp %>%
                           select(-skew, -kurtosis) %>%
                           mutate_at(vars(mean:mode),
                                     list(round),
                                     digits = 0) %>%
                           mutate_at(vars(se:upperCI),
                                     list(round),
                                     digits = 1) %>%
                           rename(estimate = mean) %>%
                           select(-median, -mode),
                         'Mark Group Site Escapement' = mark_grps %>%
                           select(-skew, -kurtosis) %>%
                           mutate_at(vars(mean:mode),
                                     list(round),
                                     digits = 0) %>%
                           mutate_at(vars(se:upperCI),
                                     list(round),
                                     digits = 1) %>%
                           rename(estimate = mean) %>%
                           select(-median, -mode)),
                    mrk_bio_list)

  #-----------------------------------------------------------------
  # write results to an Excel file
  save_list = c(list('Population Escapement' = pop_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate_at(vars(mean:mode),
                                 list(round),
                                 digits = 0) %>%
                       mutate_at(vars(se:upperCI),
                                 list(round),
                                 digits = 1) %>%
                       rename(estimate = mean) %>%
                       select(-median, -mode),
                     'All Escapement' = escape_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate_at(vars(mean:mode),
                                 list(round),
                                 digits = 0) %>%
                       mutate_at(vars(se:upperCI),
                                 list(round),
                                 digits = 1) %>%
                       rename(estimate = mean) %>%
                       select(-median, -mode),
                     'Detection' = detect_summ %>%
                       mutate_at(vars(-Node, -n_tags),
                                 list(round),
                                 digits = 3) %>%
                       rename(estimate = mean) %>%
                       select(-median, -mode),
                     'Tag Summary' = tag_summ,
                     'Biological Summary' = fullSumm),
                bio_list,
                mark_grp_list)

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
# }

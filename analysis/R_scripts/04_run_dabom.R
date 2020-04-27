# Author: Kevin See
# Purpose: prep and run DABOM
# Created: 4/1/20
# Last Modified: 4/1/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)
# library(jagsUI)
library(rjags)
library(magrittr)
library(lubridate)

spp = "Chinook"
#-----------------------------------------------------------------
# load configuration and site_df data

#-----------------------------------------------------------------
# Load required DABOM data
#-----------------------------------------------------------------
# set year
yr = 2019

load(paste0('analysis/data/derived_data/PITcleanr/TUM_', spp, '_', yr, '.rda'))

proc_ch <- proc_list$ProcCapHist %>%
  mutate(UserProcStatus = if_else(UserProcStatus == '',
                                  AutoProcStatus,
                                  UserProcStatus)) %>%
  filter(UserProcStatus)

# proc_ch %<>%
#   # filter out ghost tag detections after Sept 30
#   filter(ObsDate > lubridate::ymd(paste(yr, '0930'))) %>%
#   select(TagID) %>%
#   distinct() %>%
#   left_join(proc_ch) %>%
#   select(TagID, ObsDate, SiteID, Node, AutoProcStatus, UserProcStatus)


# # filter out ghost tag detections after Sept 30
# proc_ch %>%
#   filter(ObsDate > lubridate::ymd(paste(yr, '0930'))) %>%
#   xtabs(~ UserProcStatus, .)

# file path to the default and initial model
basic_modNm = 'analysis/model_files/TUM_DABOM.txt'

writeDABOM_TUM(file_name = basic_modNm)

#------------------------------------------------------------------------------
# Alter default model code for species and year of
# interest; sets prior for some detection node efficiencies at 0 or 100%
# based on actual tag detection data; 0% if no tags were seen
#------------------------------------------------------------------------------

# filepath for specific JAGS model code for species and year
mod_path = paste0('analysis/model_files/TUM_', spp, '_', yr, '.txt')

# writes species and year specific jags code
fixNoFishNodes(basic_modNm,
               mod_path,
               proc_ch,
               proc_list$NodeOrder)

#------------------------------------------------------------------------------
# Create capture history matrices for each main branch to be used in
# the JAGS data list
#------------------------------------------------------------------------------
# turn into wide version of capture histories, and add origin
# add fish origin from biological data
bio_df = read_rds('analysis/data/derived_data/Bio_Data_2008_2019.rds') %>%
  filter(Year == yr)

dabom_df = createDABOMcapHist(proc_ch,
                              proc_list$NodeOrder,
                              split_matrices = F) %>%
  # add origin information
  left_join(bio_df %>%
              select(TagID, Origin)) %>%
  select(TagID, Origin, everything())

dabom_list = createDABOMcapHist(proc_ch,
                                proc_list$NodeOrder,
                                split_matrices = T)

dabom_list$fishOrigin = dabom_df %>%
  mutate(Origin = recode(Origin,
                         'W' = 1,
                         'H' = 2)) %>%
  pull(Origin)

table(is.na(dabom_list$fishOrigin))

# Creates a function to spit out initial values for MCMC chains
n_branch_list = setBranchNums(parent_child)
# n_branch_list$n_pops_TUM = 8
init_fnc = setInitialValues_TUM(dabom_list,
                                n_branch_list)

# Create all the input data for the JAGS model
jags_data = createJAGSinputs_TUM(dabom_list,
                                 n_branch_list)

#------------------------------------------------------------------------------
# Tell JAGS which parameters in the model that it should save.
# the fnc is hard coded and needs to be updated if there are changes!
#------------------------------------------------------------------------------

jags_params = setSavedParams(model_file = mod_path)

#------------------------------------------------------------------------------
# Run the model

# Recommended MCMC parameters are:
# * `n.chains`: 4
# * `n.iter`: 5,000
# * `n.burnin`: 2,500
# * `n.thin`: 10
# ( 4*(5000-2500) ) / 10 = 1000 samples

set.seed(12)
# dabom_mod <- jags.basic(data = jags_data,
#                         inits = init_fnc,
#                         parameters.to.save = jags_params,
#                         model.file = mod_path,
#                         n.chains = 4,
#                         n.iter = 5000,
#                         n.burnin = 2500,
#                         n.thin = 10,
#                         # n.chains = 1,
#                         # n.iter = 2,
#                         # n.burnin = 1,
#                         # parallel = T,
#                         DIC = T)

set.seed(12)
jags = jags.model(mod_path,
                  data = jags_data,
                  inits = init_fnc,
                  n.chains = 4,
                  n.adapt = 2500)

dabom_mod = coda.samples(jags,
                         jags_params,
                         n.iter = 2500,
                         thin = 10)


# save some stuff
proc_list[["proc_ch"]] <- proc_ch

save(dabom_mod, dabom_list, proc_list, parent_child,
     file = paste0("analysis/data/derived_data/model_fits/TUM_DABOM_", spp, '_', yr,'.rda'))



#------------------------------------------------------------------------------
# diagnostics
#------------------------------------------------------------------------------
# load model run
load(paste0("analysis/data/derived_data/model_fits/PRO_DABOM_", spp, '_', yr,'.rda'))

# using mcmcr package
library(mcmcr)

# pull out mcmc.list object
my_mod = dabom_mod

#---------------------------------------
# using mcmcr
anyNA(my_mod)
my_mcmcr = as.mcmcr(my_mod)

# get Rhat statistics for all parameters
rhat_df = rhat(my_mcmcr,
               by = 'parameter',
               as_df = T) %>%
  mutate(type = if_else(grepl('_p$', parameter),
                        'Detection',
                        if_else(grepl('^p_pop', parameter) |
                                  grepl('^phi', parameter),
                                'Movement',
                                'Other')))

# plot histogram of Rhat statistics
rhat_df %>%
  ggplot(aes(x = rhat)) +
  geom_histogram(fill = 'blue',
                 bins = 40) +
  facet_wrap(~ type,
             scales = 'free')

# which parameters have converged and which haven't?
convg_df = converged(my_mcmcr,
                     by = 'parameter',
                     as_df = T)

janitor::tabyl(convg_df,
               converged)

# look at parameters that have not converged
convg_df %>%
  filter(!converged) %>%
  left_join(rhat_df)

#---------------------------------------
# using postpack
library(postpack)

# what parameters were tracked?
get_p(my_mod,
      type = 'base')

# some summary statistics
post_summ(my_mod,
          '_p$') %>%
  t() %>%
  as_tibble(rownames = 'param')

param_chk = c('p_pop_TUM')

diag_plots(post = my_mod,
           p = param_chk,
           save = T,
           file = 'outgoing/LWC_diagnostics.pdf')

# calculate Brooks-Gelman-Rubin Potential Scale Reduction Factor (Rhat)
# if ratio is close to 1, the chains have converged to the same distribution
# <1.10 is generally considered converged
post_summ(my_mod,
          # '_p$',
          get_p(my_mod,
                type = 'base'),
          ess = T, # effective sample size
          Rhat = T)[c("Rhat", "ess"),] %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(!is.na(Rhat)) %>%
  arrange(ess)

# find and remove params where Rhat == "NaN"
all_params = get_p(my_mod,
                   type = 'base')

post_summ_nas = post_summ(my_mod,
                          # '_p$',
                          all_params[-grep('deviance', all_params)],
                          ess = T, # effective sample size
                          Rhat = T)[c("Rhat", "ess"),] %>%
                  t() %>%
  as.data.frame() %>%
  as_tibble(rownames = 'param') %>%
  filter(Rhat == "NaN") %>%
  pull(param)

param_chk = get_p(my_mod, type = 'base')[grep('_p$', get_p(my_mod, type = 'base'))]
param_chk = param_chk[!param_chk %in% post_summ_nas]

# diagnostic plots for remaining params
diag_plots(post = my_mod,
           p = param_chk,
           ext_device = T)

# save plots
diag_plots(post = my_mod,
           p = param_chk,
           save = T,
           file = 'outgoing/figures/DABOM_trace_plots.pdf')

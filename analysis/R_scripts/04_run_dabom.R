# Author: Kevin See
# Purpose: prep and run DABOM
# Created: 4/1/20
# Last Modified: 5/28/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)
library(jagsUI)
# library(rjags)
library(magrittr)
library(lubridate)

#-----------------------------------------------------------------
# load configuration and site_df data

#-----------------------------------------------------------------
# Load required DABOM data
#-----------------------------------------------------------------
# set year
yr = 2012

load(paste0('analysis/data/derived_data/PITcleanr/UC_Steelhead_', yr, '.rda'))

proc_ch <- proc_list$ProcCapHist %>%
  mutate_at(vars(UserProcStatus),
            list(as.logical)) %>%
  mutate(UserProcStatus = if_else(is.na(UserProcStatus),
                                  AutoProcStatus,
                                  UserProcStatus)) %>%
  filter(UserProcStatus)

# file path to the default and initial model
basic_modNm = 'analysis/model_files/PRD_DABOM.txt'

writeDABOM_PRA(file_name = basic_modNm)

#------------------------------------------------------------------------------
# Alter default model code for species and year of
# interest; sets prior for some detection node efficiencies at 0 or 100%
# based on actual tag detection data; 0% if no tags were seen
#------------------------------------------------------------------------------

# filepath for specific JAGS model code for species and year
mod_path = paste0('analysis/model_files/PRD_Steelhead_', yr, '.txt')

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
file_nms = list.files('analysis/data/derived_data')
bio_nm = file_nms[grepl('Bio', file_nms) & grepl('.rds$', file_nms)]

bio_df = read_rds(paste0('analysis/data/derived_data/', bio_nm)) %>%
  filter(Year == yr)

dabom_df = createDABOMcapHist(proc_ch %>%
                                filter(TagID %in% bio_df$TagID),
                              proc_list$NodeOrder,
                              split_matrices = F) %>%
  # add origin information
  left_join(bio_df %>%
              select(TagID, Origin) %>%
              distinct()) %>%
  select(TagID, Origin, everything())

dabom_list = createDABOMcapHist(proc_ch %>%
                                  filter(TagID %in% bio_df$TagID),
                                proc_list$NodeOrder,
                                split_matrices = T)

dabom_list$fishOrigin = dabom_df %>%
  mutate(Origin = recode(Origin,
                         'W' = 1,
                         'H' = 2)) %>%
  pull(Origin)

table(!is.na(dabom_list$fishOrigin))

# Creates a function to spit out initial values for MCMC chains
# n_branch_list = setBranchNums(parent_child)
n_branch_list = setBranchNums_PRA()

init_fnc = setInitialValues_PRA(dabom_list)

# Create all the input data for the JAGS model
jags_data = createJAGSinputs_PRA(dabom_list)

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
dabom_mod <- jags.basic(data = jags_data,
                        inits = init_fnc,
                        parameters.to.save = jags_params,
                        model.file = mod_path,
                        n.chains = 4,
                        n.iter = 5000,
                        n.burnin = 2500,
                        n.thin = 10,
                        # n.chains = 1,
                        # n.iter = 2,
                        # n.burnin = 1,
                        # parallel = T,
                        DIC = T)

# set.seed(12)
# # adaptation phase
# jags = jags.model(mod_path,
#                   data = jags_data,
#                   inits = init_fnc,
#                   n.chains = 4,
#                   n.adapt = 1000)
# # burn-in
# update(jags,
#        n.iter = 2500)
# # posterior samples
# dabom_mod = coda.samples(jags,
#                          jags_params,
#                          n.iter = 2500,
#                          thin = 10)

# save some stuff
proc_list[["proc_ch"]] <- proc_ch

save(dabom_mod, dabom_list, proc_list, parent_child, bio_df,
     file = paste0("analysis/data/derived_data/model_fits/PRA_Steelhead_", yr,'_DABOM.rda'))



#------------------------------------------------------------------------------
# diagnostics
#------------------------------------------------------------------------------
# load model run
load(paste0("analysis/data/derived_data/model_fits/PRA_Steelhead_", yr,'_DABOM.rda'))

# using mcmcr package
library(mcmcr)

# pull out mcmc.list object
my_mod = dabom_mod

#---------------------------------------
# using mcmcr
anyNA(dabom_mod)
my_mcmcr = as.mcmcr(dabom_mod)

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

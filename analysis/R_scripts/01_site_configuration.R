# Author: Kevin See
# Purpose: Develop configuration file for DABOM
# Created: 4/1/20
# Last Modified: 6/7/21
# Notes:

# install some needed packages
install.packages(c("tidyverse",
                   "devtools",
                   "here",
                   "sf",
                   "magritter",
                   "readxl",
                   "writexl",
                   "janitor",
                   "rjags",
                   "msm",
                   "moments",
                   "coda"))

devtools::install_github("BiomarkABS/STADEM")
devtools::install_github("BiomarkABS/PITcleanr")
devtools::install_github("BiomarkABS/DABOM")

#------------------------------,-----------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(sf)
library(here)

#-----------------------------------------------------------------
# set starting point
root_site = "PRA"

# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  # manually add site for Colockum Creek (not in PTAGIS)
  bind_rows(tibble(site_code = 'CLK',
                   config_id = 100,
                   antenna_id = 'A1',
                   node = 'CLK',
                   # making these up
                   start_date = as.POSIXct(lubridate::ymd('20100101')),
                   site_type = 'INT',
                   site_name = 'Colockum Creek',
                   antenna_group = 'Single Colockum Ck',
                   site_description = 'Temporary single antenna.',
                   site_type_name = 'Instream Remote Detection System',
                   rkm = '740.001',
                   rkm_total = 741,
                   # this puts CLK upstream of RIA
                   latitude = 47.3707357269787,
                   longitude = -120.25617371760839)) %>%
  # this puts CLK on Colockum Creek, but between PRA and RIA
  # latitude = 47.29722788926544,
  # longitude = -120.10577913008702)) %>%
  filter(!(site_code == 'WAN' & site_type == 'MRR'),
         !(site_code == 'TMF' & site_type == 'MRR'),
         !(site_code == 'PRO' & site_type == 'MRR')) %>%
  mutate(node = if_else(site_code %in% c('RIA', 'RRF', 'WEA', 'PRV'),
                        site_code,
                        node)) %>%
  mutate(node = if_else(site_code == 'PRDLD1',
                        'PRA',
                        node)) %>%
  mutate(node = if_else(node == "LWE",
                        'LWEB0',
                        node),
         node = if_else(site_code %in% c('TUF', 'TUMFBY', 'TUM'),
                        'TUM',
                        node),
         node = if_else(site_code == 'LNF' & antenna_id %in% c('01', '02'),
                        'LNFA0',
                        node),
         node = if_else(site_code == 'LNF' & antenna_id %in% c('03', '04'),
                        'LNFB0',
                        node),
         node = if_else(site_code == 'LEAV',
                        'LNFA0',
                        node),
         node = if_else(site_code == 'ICL' & config_id == 100,
                        'ICLB0',
                        node),
         node = if_else(site_code == 'CHIWAC',
                        'CHWA0',
                        node),
         node = if_else(site_code == 'CHIWAR',
                        'CHLA0',
                        node),
         node = if_else(site_code == 'CHIKAC',
                        'CHUA0',
                        node),
         node = if_else(site_code == 'WHITER',
                        'WTLA0',
                        node),
         node = if_else(site_code == 'LWENAT',
                        'LWNA0',
                        node),
         node = if_else(site_code == 'NASONC',
                        'NALA0',
                        node),
         # any fish seen at Dryden dam should also be seen at LWE
         node = if_else(site_code == 'DRY',
                        'LWEA0',
                        node),
         # any fish seen at Chiwawa acclimation pond gets moved to CHL
         node = if_else(site_code == 'CHP',
                        'CHLA0',
                        node),
         node = if_else(site_code == 'EBO',
                        'RRF',
                        node),
         node = if_else(site_code == 'RRJ',
                        'RRF',
                        node),
         node = if_else(site_code == 'EHL' & config_id == 100 & antenna_id == '02',
                        'EHLB0',
                        node),
         node = if_else(site_code == 'EHL' & config_id == 100 & antenna_id == '01',
                        'EHLA0',
                        node),
         node = if_else(site_code == 'EHL' & config_id == 110 & antenna_id == '03',
                        'EHLB0',
                        node),
         node = if_else(site_code == 'EHL' & config_id == 110 & antenna_id %in% c('01', '02'),
                        'EHLA0',
                        node),
         # combine a couple sites in the Entiat
         node = if_else(site_code %in% c("ENS", "ENM"),
                        "ENAA0",
                        node),
         node = if_else(site_code == "WEH" & antenna_id == "A2",
                        "WEHB0",
                        node),
         node = if_else(site_code == "WEH" & antenna_id != "A2",
                        "WEHA0",
                        node),
         node = if_else(node == "LMR",
                        'LMRB0',
                        node),
         node = if_else(site_code == 'LBC' & config_id == 100,
                        'LBCB0',
                        node),
         node = if_else(site_code == 'MRC',
                        'MRCB0',
                        node),
         node = if_else(site_code %in% c('SSC', '18N', 'MHB', 'M3R', 'MWF'),
                        'MRCA0',
                        node),
         node = if_else(site_code == 'MSH' & antenna_id %in% c('02', '03'),
                        'MSHB0',
                        node),
         node = if_else(site_code == 'MSH' & antenna_id %in% c('01'),
                        'MSHA0',
                        node),
         node = if_else(site_code == 'MSH' & antenna_id == '00',
                        'METHB0',
                        node),
         node = if_else(site_code == 'METH',
                        'METHA0',
                        node),
         node = if_else(site_code == 'LLC' & config_id == 100,
                        if_else(antenna_id == 'D3',
                                'LLCB0',
                                'LLCA0'),
                        node),
         node = if_else(node == "SCP",
                        'SCPB0',
                        node),
         # node = if_else(node == "OMK",
         #               'OMKB0',
         #               node),
         # node = if_else(site_code %in% c('OBF', 'OMF'),
         #               'OMKA0',
         #               node),
         node = if_else(site_code == "OMF",
                        "OBF",
                        node),
         node = if_else(site_code == 'ZSL',
                        if_else(grepl('Weir 3', antenna_group, ignore.case = T),
                                'ZSLB0',
                                'ZSLA0'),
                        node),
         node = if_else(site_code == 'SA1' & config_id == 110,
                        'SA1B0',
                        node),
         node = if_else(site_code == 'OKC' & config_id == 100,
                        'OKCB0',
                        node),
         node = if_else(site_code == 'OMK' & config_id == 100,
                        'OMKB0',
                        node),
         # combine some sites above OKV into the upstream array at OKV
         node = if_else(site_code %in% c("OKS", "OKW"),
                        "OKVA0",
                        node),
         node = if_else(site_code == 'RCT' & config_id == 100,
                        'RCTB0',
                        node),
         node = if_else(site_code == 'BPC' & config_id == 100,
                        if_else(antenna_id %in% c('C3'),
                                'BPCB0',
                                'BPCA0'),
                        node),
         node = if_else(site_code == 'PRH' & antenna_id %in% c('F1', 'F2', 'F3', 'F4'),
                        'PRHB0',
                        node),
         node = if_else((site_code == 'PRH' & antenna_id %in% c('F5', 'F6', '01', '02')) | site_code %in% c('DDM', 'DM', 'UM', 'UUM', 'UP'),
                        'PRHA0',
                        node),
         node = if_else(site_code == 'PRO' & site_type == 'INT',
                        'PROB0',
                        node),
         # grab all sites upstream of Prosser dam, and assign them to PROA0
         node = if_else(site_code != "PRO" &
                          as.integer(stringr::str_split(rkm, '\\.', simplify = T)[,1]) == 539 &
                          as.integer(stringr::str_split(rkm, '\\.', simplify = T)[,2]) > 76,
                        "PROA0",
                        node),
         node = if_else(site_code == 'ICH',
                        'ICHB0',
                        node),
         node = if_else(grepl('522\\.', rkm) & rkm_total > 538,
                        'ICHA0',
                        node),
         node = if_else(site_code == 'MDR',
                        'MDRB0',
                        node),
         node = if_else(site_code %in% c('LWD', 'BGM', 'NBA', 'MCD'),
                        'MDRA0',
                        node),
         node = if_else(site_code == 'HST',
                        'HSTB0',
                        node),
         node = if_else(site_code %in% c('BBT', 'COP', 'PAT'),
                        'HSTA0',
                        node),
         node = if_else(site_code == 'JD1',
                        'JD1B0',
                        node),
         node = if_else(site_code %in% c('30M', 'BR0', 'JDM', 'SJ1', 'SJ2', 'MJ1'),
                        'JD1A0',
                        node),
         node = if_else(site_code != 'JD1' & as.integer(stringr::str_split(rkm, '\\.', simplify = T)[,1]) < 351,
                        'JDA',
                        node)) %>%
  distinct() %>%
  # correct a couple rkm values
  mutate(rkm = if_else(site_code == 'SA1',
                       '858.041.003',
                       rkm),
         rkm_total = if_else(site_code == 'SA1',
                             902,
                             rkm_total)) %>%
  mutate(rkm = if_else(site_code == 'TON',
                       '858.133.001',
                       rkm),
         rkm_total = if_else(site_code == 'TON',
                             992,
                             rkm_total)) %>%
  mutate(rkm = if_else(grepl('WEH', node),
                       '829.001',
                       rkm),
         rkm_total = if_else(grepl('WEH', node),
                             830,
                             rkm_total)) %>%
  mutate(rkm = if_else(site_code == "MSH",
                       '843.082',
                       rkm),
         rkm_total = if_else(site_code == "MSH",
                             925,
                             rkm_total),
         rkm = if_else(site_code == "METH",
                       '843.083',
                       rkm),
         rkm_total = if_else(site_code == "METH",
                             926,
                             rkm_total))

# Node network for DABOM

# get spatial object of sites used in model
sites_sf = writeOldNetworks()$PriestRapids %>%
  mutate(across(c(SiteID, Step3),
                recode,
                "BelowJD1" = "JDA"),
         path = str_replace(path, "BelowJD1", "JDA")) %>%
  rename(site_code = SiteID) %>%
  left_join(configuration) %>%
  group_by(site_code) %>%
  filter(config_id == max(config_id)) %>%
  ungroup() %>%
  select(site_code,
         site_name,
         site_type = site_type_name,
         type = site_type,
         rkm,
         site_description = site_description,
         latitude, longitude) %>%
  distinct() %>%
  filter(!is.na(latitude)) %>%
  st_as_sf(coords = c("longitude",
                      "latitude"),
           crs = 4326) %>%
  st_transform(crs = 5070)

#-----------------------------------------------------------------
# download the NHDPlus v2 flowlines
# do you want flowlines downstream of root site? Set to TRUE if you have downstream sites
dwn_flw = T
nhd_list = queryFlowlines(sites_sf = sites_sf,
                          root_site_code = root_site,
                          min_strm_order = 2,
                          dwnstrm_sites = dwn_flw,
                          dwn_min_stream_order_diff = 4)

# compile the upstream and downstream flowlines
flowlines = nhd_list$flowlines
if(dwn_flw) {
  flowlines %<>%
    rbind(nhd_list$dwn_flowlines)
}

# upstream extent of study area (cut off areas further upstream)
upstrm_loc = "Chief Joseph Dam"

library(ggmap)

upstrm_comid = ggmap::geocode(upstrm_loc, output = "latlon") %>%
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326) %>%
  nhdplusTools::discover_nhdplus_id()

nhd_upstrm_lst = nhdplusTools::plot_nhdplus(outlets = list(upstrm_comid),
                                            streamorder = min(nhd_list$flowlines$StreamOrde),
                                            actually_plot = F)

flowlines %<>%
  anti_join(nhd_upstrm_lst$flowline %>%
              st_drop_geometry() %>%
              select(Hydroseq))
#-----------------------------------------------------------------
# plot the flowlines and the sites
ggplot() +
  geom_sf(data = flowlines,
          aes(color = as.factor(StreamOrde),
              size = StreamOrde)) +
  scale_color_viridis_d(direction = -1,
                        option = "D",
                        end = 0.8) +
  scale_size_continuous(range = c(0.2, 1.2),
                        guide = 'none') +
  # geom_sf(data = nhd_list$basin,
  #         fill = NA,
  #         lwd = 2) +
  # this cuts out parts of the basin upstream of upstrm_loc
  geom_sf(data = flowlines %>%
            filter(!Hydroseq %in% nhd_list$dwn_flowlines$Hydroseq) %>%
            summarise(bndry = 'basin') %>%
            select(bndry) %>%
            st_convex_hull(),
          fill = NA,
          lwd = 2) +
  geom_sf(data = sites_sf,
          size = 4,
          color = "black") +
  geom_sf_label(data = sites_sf,
                aes(label = site_code)) +
  geom_sf_label(data = sites_sf %>%
                  filter(site_code == root_site),
                aes(label = site_code),
                color = "red") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(color = "Stream\nOrder")


#-----------------------------------------------------------------
# build parent child table
parent_child = sites_sf %>%
  filter(! site_code %in% c("EBO")) %>%
  # filter(! site_code %in% c("MC1", "MC2", "MCJ",
  #                           "FDD",
  #                           "MWC",
  #                           "YHC",
  #                           "HN1", "HN3",
  #                           "HSM", "HSU",
  #                           "LOR",
  #                           "EWC",
  #                           "LBT",
  #                           "LTP",
  #                           "BBP", "BCP",
  #                           "MRT", "HSL",
  #                           "ENM", 'ENS', 'TY4', '3D4')) %>%
  buildParentChild(flowlines,
                   # rm_na_parent = T,
                   add_rkm = F) %>%
  editParentChild(fix_list = list(c("JDA", 'ICH', "PRA"),
                                  c("JDA", 'RSH', "PRA"),
                                  c("JDA", 'JD1', "PRA"),
                                  c("JDA", 'PRO', "PRA"),
                                  c("JDA", 'TMF', "PRA"),
                                  c("JDA", 'PRV', "PRA"),
                                  c(NA, "JDA", 'PRA'),
                                  c("RSH", 'PRH', 'PRA'),
                                  c("ICL", 'TUM', "LWE"),
                                  c("LNF", 'ICM', "ICL"),
                                  c(NA, "LNF", "ICL"),
                                  c("RIA", "WEA", 'RRF'),
                                  c("RIA", "WEH", 'RRF'),
                                  c("RIA", "ENL", "RRF"),
                                  c("WEH", "LMR", "WEA"),
                                  c("WEH", "OKL", "WEA"),
                                  c("WEH", "FST", 'WEA'),
                                  c("EHL", 'ENA', 'ENL'),
                                  c("EHL", 'MAD', 'ENL'),
                                  c("METH", "MRW", "MRC"),
                                  c("SCP", "METH", "MSH"),
                                  c("SCP", 'MSH', 'MRC'),
                                  c("WHS", "OKC", "ZSL"),
                                  c("WHS", "ZSL", "OKL"),
                                  c("ZSL", 'OKV', 'OKC'),
                                  c("JOH", 'WHS', 'OKL'),
                                  c("JOH", 'BPC', 'OKL'),
                                  c("JOH", 'ANT', 'OKL'),
                                  c("JOH", 'TNK', 'OKL'),
                                  c("JOH", 'AEN', 'OKL')),
                  switch_parent_child = list(c("RSH", "PRA"))) %>%
  filter(!parent %in% c("WEH", "PRH"))

# add RKMs from configuration file (since we had to fix at least one from PTAGIS)
parent_child %<>%
  left_join(configuration %>%
              select(parent = site_code,
                     parent_rkm = rkm) %>%
              distinct(),
            by = "parent") %>%
  left_join(configuration %>%
              select(child = site_code,
                     child_rkm = rkm) %>%
              distinct(),
            by = "child") %>%
  distinct()


sites_df = writeOldNetworks()$PriestRapids %>%
  mutate(across(c(SiteID, Step3),
                recode,
                "BelowJD1" = "JDA"),
         path = str_replace(path, "BelowJD1", "JDA")) %>%
  rename(site_code = SiteID)

ques_locs = sites_df %>%
  # filter(grepl('Wenatchee', path)) %>%
  # filter(grepl("Entiat", path)) %>%
  # filter(grepl('Methow', path)) %>%
  # filter(grepl("Okanogan", path)) %>%
  filter(nchar(Step3) == 3,
         nchar(Step4) < 6,
         nchar(Step5) < 6,
         Step2 != "BelowPriest") %>%
  pull(site_code)

parent_child %>%
  filter(parent %in% ques_locs |
           child %in% ques_locs) %>%
  buildPaths() %>%
  left_join(sites_df %>%
              select(end_loc = site_code,
                     org_path = path))

#-----------------------------------------------------------------
# Save file.
save(configuration,
     sites_sf,
     flowlines,
     parent_child,
     file = here('analysis/data/derived_data/site_config.rda'))


#-----------------------------------------------------------------
# Build network diagram
# simple
pc_graph = plotNodes(parent_child,
                     layout = "tree")

pc_nodes_graph = parent_child %>%
  addParentChildNodes(configuration) %>%
  plotNodes()

# control more settings
node_order = buildNodeOrder(parent_child) %>%
  separate(col = path,
           into = paste("step", 1:max(.$node_order), sep = "_"),
           remove = F) %>%
  mutate(branch_nm = if_else(node == "PRA",
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
                                                                     "Downstream",
                                                                     "Other"))))))) %>%
  select(-starts_with("step"))

nodes = buildNodeGraph(parent_child) %>%
  as_tibble() %>%
  left_join(node_order %>%
              select(label = node,
                     branch_nm))
edges = parent_child %>%
  left_join(nodes, by = c('parent' = 'label')) %>%
  rename(from = index) %>%
  left_join(nodes, by = c('child' = 'label')) %>%
  rename(to = index) %>%
  select(from, to)

library(ggraph)
node_graph = tidygraph::tbl_graph(nodes = nodes,
                                  edges = edges)

node_p = node_graph %>%
  ggraph(layout = "tree") +
  # ggraph(layout = "partition") +
  # ggraph(layout = "kk") +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm'),
                               type = "closed"),
                 end_cap = circle(4, 'mm')) +
  geom_node_point(size = 7,
                  aes(color = branch_nm)) +
  theme_graph(base_family = 'Times') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette = "Set1",
                     na.value = "black") +
  geom_node_label(aes(label = label),
                  size = 2,
                  label.padding = unit(0.1, 'lines'),
                  label.size = 0.1)

node_p

# save as pdf
library(here)
ggsave(here("analysis/figures/PriestRapids_DABOM_sites.pdf"),
       node_p,
       width = 9,
       height = 6)

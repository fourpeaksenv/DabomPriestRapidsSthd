# Author: Kevin See
# Purpose: Develop configuration file for DABOM
# Created: 4/27/20
# Last Modified: 12/17/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(magrittr)

#-----------------------------------------------------------------
# dataframe of sites for PRD DABOM model with some indication of network structure
site_df = writePRDNodeNetwork()

# build configuration table (requires internet connection)
org_config = buildConfig()

# manually add site for Colockum Creek (not in PTAGIS)
org_config = org_config %>%
  bind_rows(tibble(SiteID = 'CLK',
                   ConfigID = 100,
                   AntennaID = 'A1',
                   Node = 'CLK',
                   ValidNode = T,
                   # making these up
                   StartDate = as.POSIXct(lubridate::ymd('20100101')),
                   SiteType = 'INT',
                   SiteName = 'Colockum Creek',
                   AntennaGroup = 'Single Colockum Ck',
                   SiteDescription = 'Tempoary single antenna.',
                   SiteTypeName = 'Instream Remote Detection System',
                   RKM = '740.001',
                   RKMTotal = 741))

# customize some nodes based on DABOM framework
configuration = org_config %>%
  filter(!(SiteID == 'WAN' & SiteType == 'MRR'),
         !(SiteID == 'TMF' & SiteType == 'MRR'),
         !(SiteID == 'PRO' & SiteType == 'MRR')) %>%
  mutate(Node = if_else(SiteID %in% c('RIA', 'RRF', 'WEA', 'PRV'),
                       SiteID,
                       Node)) %>%
  mutate(Node = if_else(SiteID == 'PRDLD1',
                       'PRA',
                       Node)) %>%
  mutate(Node = if_else(Node == "LWE",
                       'LWEB0',
                       Node),
         Node = if_else(SiteID %in% c('TUF', 'TUMFBY', 'TUM'),
                       'TUM',
                       Node),
         Node = if_else(SiteID == 'LNF' & AntennaID %in% c('01', '02'),
                       'LNFA0',
                       Node),
         Node = if_else(SiteID == 'LNF' & AntennaID %in% c('03', '04'),
                       'LNFB0',
                       Node),
         Node = if_else(SiteID == 'LEAV',
                       'LNFA0',
                       Node),
         Node = if_else(SiteID == 'ICL' & ConfigID == 100,
                       'ICLB0',
                       Node),
         Node = if_else(SiteID == 'CHIWAC',
                       'CHWA0',
                       Node),
         Node = if_else(SiteID == 'CHIWAR',
                       'CHLA0',
                       Node),
         Node = if_else(SiteID == 'CHIKAC',
                       'CHUA0',
                       Node),
         Node = if_else(SiteID == 'WHITER',
                       'WTLA0',
                       Node),
         Node = if_else(SiteID == 'LWENAT',
                       'LWNA0',
                       Node),
         Node = if_else(SiteID == 'NASONC',
                       'NALA0',
                       Node),
         # any fish seen at Dryden dam should also be seen at LWE
         Node = if_else(SiteID == 'DRY',
                       'LWEA0',
                       Node),
         # any fish seen at Chiwawa acclimation pond gets moved to CHL
         Node = if_else(SiteID == 'CHP',
                       'CHLA0',
                       Node),
         Node = if_else(SiteID == 'EBO',
                       'RRF',
                       Node),
         Node = if_else(SiteID == 'EHL' & ConfigID == 100 & AntennaID == '02',
                       'EHLB0',
                       Node),
         Node = if_else(SiteID == 'EHL' & ConfigID == 100 & AntennaID == '01',
                       'EHLA0',
                       Node),
         Node = if_else(SiteID == 'EHL' & ConfigID == 110 & AntennaID == '03',
                       'EHLB0',
                       Node),
         Node = if_else(SiteID == 'EHL' & ConfigID == 110 & AntennaID %in% c('01', '02'),
                       'EHLA0',
                       Node),
         Node = if_else(SiteID == "WEH" & AntennaID == "A2",
                        "WEHB0",
                        Node),
         Node = if_else(SiteID == "WEH" & AntennaID != "A2",
                        "WEHA0",
                        Node),
         Node = if_else(Node == "LMR",
                       'LMRB0',
                       Node),
         Node = if_else(SiteID == 'LBC' & ConfigID == 100,
                       'LBCB0',
                       Node),
         Node = if_else(SiteID == 'MRC',
                       'MRCB0',
                       Node),
         Node = if_else(SiteID %in% c('SSC', '18N', 'MHB', 'M3R', 'MWF'),
                       'MRCA0',
                       Node),
         Node = if_else(SiteID == 'MSH' & AntennaID %in% c('02', '03'),
                       'MSHB0',
                       Node),
         Node = if_else(SiteID == 'MSH' & AntennaID %in% c('01'),
                       'MSHA0',
                       Node),
         Node = if_else(SiteID == 'MSH' & AntennaID == '00',
                       'METHB0',
                       Node),
         Node = if_else(SiteID == 'METH',
                       'METHA0',
                       Node),
         Node = if_else(SiteID == 'LLC' & ConfigID == 100,
                       if_else(AntennaID == 'D3',
                              'LLCB0',
                              'LLCA0'),
                       Node),
         Node = if_else(Node == "SCP",
                       'SCPB0',
                       Node),
         Node = if_else(Node == "OMK",
                       'OMKB0',
                       Node),
         Node = if_else(SiteID %in% c('OBF', 'OMF'),
                       'OMKA0',
                       Node),
         Node = if_else(SiteID == 'ZSL',
                       if_else(grepl('Weir 3', AntennaGroup, ignore.case = T),
                              'ZSLB0',
                              'ZSLA0'),
                       Node),
         Node = if_else(SiteID == 'SA1' & ConfigID == 110,
                       'SA1B0',
                       Node),
         Node = if_else(SiteID == 'OKC' & ConfigID == 100,
                       'OKCB0',
                       Node),
         # combine some sites above OKV into the upstream array at OKV
         Node = if_else(SiteID %in% c("OKS", "OKW"),
                        "OKVA0",
                        Node),
         Node = if_else(SiteID == 'RCT' & ConfigID == 100,
                       'RCTB0',
                       Node),
         Node = if_else(SiteID == 'BPC' & ConfigID == 100,
                       if_else(AntennaID %in% c('C3'),
                              'BPCB0',
                              'BPCA0'),
                       Node),
         Node = if_else(SiteID == 'PRH' & AntennaID %in% c('F1', 'F2', 'F3', 'F4'),
                       'PRHB0',
                       Node),
         Node = if_else((SiteID == 'PRH' & AntennaID %in% c('F5', 'F6', '01', '02')) | SiteID %in% c('DDM', 'DM', 'UM', 'UUM', 'UP'),
                       'PRHA0',
                       Node),
         Node = if_else(SiteID == 'PRO' & SiteType == 'INT',
                       'PROB0',
                       Node),
         Node = if_else(SiteID %in% c('CHANDL', 'SAT', 'TOP', 'SUN', 'LNR', 'ROZ', 'LMC', 'TAN') | SiteID == 'PRO' & SiteType == 'MRR',
                       'PROA0',
                       Node),
         Node = if_else(SiteID == 'ICH',
                       'ICHB0',
                       Node),
         Node = if_else(grepl('522\\.', RKM) & RKMTotal > 538,
                       'ICHA0',
                       Node),
         Node = if_else(SiteID == 'MDR',
                       'MDRB0',
                       Node),
         Node = if_else(SiteID %in% c('LWD', 'BGM', 'NBA', 'MCD'),
                       'MDRA0',
                       Node),
         Node = if_else(SiteID == 'HST',
                       'HSTB0',
                       Node),
         Node = if_else(SiteID %in% c('BBT', 'COP', 'PAT'),
                       'HSTA0',
                       Node),
         Node = if_else(SiteID == 'JD1',
                       'JD1B0',
                       Node),
         Node = if_else(SiteID %in% c('30M', 'BR0', 'JDM', 'SJ1', 'SJ2', 'MJ1'),
                       'JD1A0',
                       Node),
         Node = if_else(SiteID != 'JD1' & as.integer(stringr::str_split(RKM, '\\.', simplify = T)[,1]) < 351,
                       'BelowJD1',
                       Node)) %>%
  distinct()

# correct a couple RKM values
configuration = configuration %>%
  mutate(RKM = if_else(SiteID == 'SA1',
                      '858.041.003',
                      RKM),
         RKMTotal = if_else(SiteID == 'SA1',
                           902,
                           RKMTotal)) %>%
  mutate(RKM = if_else(SiteID == 'TON',
                      '858.133.001',
                      RKM),
         RKMTotal = if_else(SiteID == 'TON',
                           992,
                           RKMTotal)) %>%
  mutate(RKM = if_else(grepl('WVT', Node),
                      '829.001',
                      RKM),
         RKMTotal = if_else(grepl('WVT', Node),
                           830,
                           RKMTotal))

# group sites appropriately
site_list = vector('list', 5)
names(site_list) = c('BelowPriest', 'Wenatchee', 'Entiat', 'Methow', 'Okanogan')
for(grp in names(site_list)) {
  site_list[[grp]] = site_df %>%
    filter(grepl(grp, path)) %>%
    select(SiteID) %>%
    as.matrix() %>%
    as.character()
}
site_list[['Wenatchee']] = c('RIA', site_list[['Wenatchee']])
site_list[['Entiat']] = c('RRF', 'WEH', site_list[['Entiat']])
site_list[['Methow']] = c('WEA', site_list[['Methow']])

# Save file.
save(configuration,
     site_df,
     site_list,
     file = 'analysis/data/derived_data/site_config.rda')


#-----------------------------------------------------------------
# which sites are in site_df, but not in the PTAGIS configuration file?
site_df %>%
  filter(!(SiteID %in% configuration$SiteID |
             SiteID %in% configuration$Node))

#-----------------------------------------------------------------
# Build network diagram
#-----------------------------------------------------------------
library(tidygraph)
library(ggraph)
library(netplot)

# build parent-child table
# which spawn year are we dealing with?
yr = 2020
# start date is July 1 of the previous year
start_date = paste0(yr - 1, '0701')

# build parent-child table
par_ch_node = createParentChildDf(site_df,
                                  configuration,
                                  startDate = start_date)

root_node = par_ch_node %>%
  filter(nodeOrder == 1) %>%
  pull(ParentNode)

par_ch_site = createParentChildDf(site_df,
                                  configuration %>%
                                    mutate(Node = if_else(grepl('A0$', Node) | grepl('B0$', Node),
                                                         SiteID,
                                                         Node)) %>%
                                    distinct(),
                                  startDate = start_date) %>%
  rename(ParentSite = ParentNode,
         ChildSite = ChildNode) %>%
  left_join(stack(site_list) %>%
              tbl_df() %>%
              select(Group = ind,
                     ChildSite = values) %>%
              bind_rows(tibble(Group = root_node,
                               ChildSite = root_node)) %>%
              mutate(Group = factor(Group,
                                    levels = c(root_node, names(site_list)))) %>%
              mutate(BranchNum = as.integer(Group))) %>%
  left_join(configuration %>%
              select(ChildSite = SiteID,
                     lat = Latitude,
                     long = Longitude) %>%
              distinct()) %>%
  mutate(lat = if_else(is.na(lat) & ChildSite == 'BelowJD1',
                       configuration %>%
                         filter(SiteID == 'CHINOR') %>%
                         pull(Latitude),
                       lat),
         long = if_else(is.na(long) & ChildSite == 'BelowJD1',
                        configuration %>%
                          filter(SiteID == 'CHINOR') %>%
                          pull(Longitude),
                        long))

# add nodes for black boxes
bb_nodes = par_ch_site %>%
  group_by(ParentSite) %>%
  summarise(nChild = n_distinct(ChildSite)) %>%
  filter(nChild > 1) %>%
  left_join(par_ch_site %>%
              select(-ParentSite,
                     ParentSite = ChildSite,
                     SiteType:long)) %>%
  mutate(SiteType = 'BB') %>%
  mutate(ChildSite = paste0(ParentSite, '_bb')) %>%
  select(-nChild) %>%
  distinct()

par_ch_site %<>%
  bind_rows(bb_nodes)

# build table of nodes
nodes = par_ch_site %>%
  select(ParentSite, ChildSite) %>%
  gather(type, node) %>%
  select(node) %>%
  distinct() %>%
  left_join(par_ch_site %>%
              rename(node = ChildSite) %>%
              select(-starts_with("Parent"))) %>%
  mutate(Group = as.factor(Group),
         Group = fct_relevel(Group, root_node)) %>%
  arrange(Group, RKM, nodeOrder) %>%
  mutate(index = 1:n()) %>%
  select(index, label = node, everything())

nodes %<>%
  mutate(Group = if_else(label %in% c('PRA', 'PRA_bb', 'RIA', 'RRF', 'WEA'),
                         'Mainstem',
                         as.character(Group)),
         Group = factor(Group,
                        levels = c('Mainstem',
                                   'BelowPriest',
                                   'Wenatchee',
                                   'Entiat',
                                   'Methow',
                                   'Okanogan')),
         Group = fct_drop(Group))

nodes %<>%
  left_join(par_ch_site %>%
              group_by(label = ParentSite) %>%
              summarise(nChilds = n_distinct(ChildSite)) %>%
              bind_rows(par_ch_site %>%
                          filter(!ChildSite %in% ParentSite) %>%
                          select(label = ChildSite) %>%
                          mutate(nChilds = 0)) %>%
              mutate(nodeType = if_else(nChilds == 0,
                                        'Terminal',
                                        if_else(nChilds == 1,
                                                'PassThru', 'Branch')))) %>%
  mutate(nodeType = if_else(SiteType == 'BB',
                            'BB', nodeType),
         nodeType = factor(nodeType,
                           levels = c('Branch',
                                      'PassThru',
                                      'Terminal',
                                      'BB')))


# build table of edges (connecting nodes)
edges = par_ch_site %>%
  filter(ParentSite != ChildSite) %>%
  select(from = ParentSite,
         to = ChildSite) %>%
  distinct() %>%
  mutate(edgeID = 1:n()) %>%
  gather(direction, label, -edgeID) %>%
  left_join(nodes %>%
              select(index, label)) %>%
  select(-label) %>%
  spread(direction, index) %>%
  select(-edgeID)

# one graph with all sites
myGraph = tbl_graph(nodes = nodes,
                    edges = edges)


#--------------------------------------------------
# tidygraph

# myLayout = c('kk')
# myLayout = c('dendrogram')
# myLayout = c('treemap')
myLayout = c('partition')
# myLayout = c('circlepack')
myLayout = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds',
             'randomly', 'fr', 'kk', 'drl', 'lgl')[9]
#
# c('star', 'dh', 'lgl')

set.seed(8)
allGr_p = myGraph %>%
  # ggraph(layout = 'igraph',
  #        algorithm = 'nicely') +
  # ggraph(layout = 'igraph',
  #        algorithm = 'fr') +
  ggraph(layout = myLayout) +
  geom_edge_diagonal() +
  # geom_edge_link() +
  # geom_node_point(aes(color = Group,
  #                     shape = nodeType),
  #                 size = 7) +
  geom_node_point(aes(color = Group),
                  size = 7) +
  geom_node_label(aes(label = label),
                  # repel = T,
                  size = 2,
                  label.padding = unit(0.1, 'lines'),
                  label.size = 0.1) +
  # geom_node_text(aes(label = label),
  #                 # repel = T,
  #                 size = 2) +
  scale_color_brewer(palette = 'Set1',
                     guide = 'none') +
  scale_shape_manual(values = c('Branch' = 19,
                                'PassThru' = 18,
                                'Terminal' = 15,
                                'BB' = 22),
                     guide = 'none') +
  theme_graph(base_family = 'Times') +
  theme(legend.position = 'bottom')
allGr_p

ggsave('analysis/figures/UC_SiteSchematic_2020.pdf',
       allGr_p,
       width = 15,
       height = 5,
       dpi = 600)


#--------------------------------------------------
# netplot
set.seed(8)
# l = igraph::layout_with_fr(myGraph)
# l = igraph::layout_with_kk(myGraph)
l = igraph::layout_with_dh(myGraph)
# l = igraph::layout_with_gem(myGraph)
# l = igraph::layout_nicely(myGraph)
# l = igraph::layout_on_grid(myGraph)
# l = igraph::layout_in_circle(myGraph)
# l = igraph::layout_with_graphopt(myGraph)
# l = igraph::layout_with_sugiyama(myGraph)
l = igraph::layout_as_tree(myGraph,
                           circular = T,
                           flip.y = F)

# l = igraph::layout_with_lgl(myGraph,
#                             root = 1)
#
# l = igraph::layout_with_mds(myGraph,
#                             dist = dist(st_coordinates(nodes_sf),
#                                         method = 'euclidean',
#                                         diag = T,
#                                         upper = T) %>%
#                               as.matrix(),
#                             dim = 2)

# set of colors
myColors = RColorBrewer::brewer.pal(nlevels(nodes$Group), 'Set1')
# myColors = viridis::plasma(nlevels(nodes$Group))
# myColors = gray.colors(nlevels(nodes$Group), start = 0.3, end = 0.9)
# myColors = c(rep('darkgray', 3), 'black')

# myColors = sample(myColors, length(myColors))

# myLabs = nodes$label
# myLabs[grepl('_bb$', myLabs)] = NA

prd_sites = nplot(myGraph,
                  layout = l,
                  vertex.color = myColors[nodes$Group],
                  # vertex.color = myColors[nodes$nodeType],
                  vertex.size = 1,
                  # vertex.size.range = c(0.1, 0.03),
                  # vertex.nsides = c(100, 3, 4, 4)[nodes$nodeType],
                  vertex.nsides = 100,
                  # vertex.rot = c(0, 1.55, 0.78, 0.78)[nodes$nodeType],
                  vertex.label = nodes$label,
                  vertex.label.fontsize = 5,
                  edge.curvature = pi/6,
                  zero.margins = T)
prd_sites

pdf('analysis/figures/UC_SiteSchematic_v3.pdf',
    width = 8,
    height = 6)
print(prd_sites)
dev.off()

#--------------------------------------------------
# igraph
library(igraph)
l = igraph::layout_as_tree(myGraph,
                           flip.y = F)

# set of colors
myColors = RColorBrewer::brewer.pal(nlevels(nodes$Group), 'Set1')
# myColors = viridis::plasma(nlevels(nodes$Group))
# myColors = viridis::viridis(nlevels(nodes$Group))
# myColors = gray.colors(nlevels(nodes$Group), start = 0.3, end = 0.9)
# myColors = c(rep('darkgray', 3), 'black')

# this will open up a Quartz window and let you edit the layout with your mouse
id = tkplot(myGraph,
            layout = l,
            canvas.width = 1400,
            canvas.height = 800,
            vertex.color = myColors[nodes$Group],
            # vertex.color = myColors[nodes$nodeType],
            # vertex.color = 'gray90',
            vertex.shape = c('circle', 'square')[1],
            # vertex.shape = c('circle', 'diamond', 'square', 'square')[nodes$nodeType],
            vertex.size = 16,
            # label = nodes$label,
            label.size = 1,
            label.color = 'red',
            edge.color = 'black')

# for saving
tk_postscript(tkp.id = id)

tk_close(tkp.id = id)
tk_off()


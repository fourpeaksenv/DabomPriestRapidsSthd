# Author: Kevin See
# Purpose: Develop configuration file for DABOM
# Created: 4/1/20
# Last Modified: 4/1/20
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(magrittr)

#-----------------------------------------------------------------
# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  mutate(Node = ifelse(SiteID %in% c('LNF', 'LEAV'),
                       'LNF',
                       Node),
         Node = ifelse(SiteID %in% c('TUF', 'TUMFBY', 'TUM'),
                       'TUM',
                       Node),
         Node = ifelse(SiteID == 'CHIWAC',
                       'CHWA0',
                       Node),
         Node = ifelse(SiteID == 'CHIWAR',
                       'CHLA0',
                       Node),
         Node = ifelse(SiteID == 'CHIW',
                       'CHLA0',
                       Node),
         Node = ifelse(SiteID == 'CHIKAC',
                       'CHUA0',
                       Node),
         Node = ifelse(SiteID == 'WHITER',
                       'WTLA0',
                       Node),
         Node = ifelse(SiteID == 'LWENAT',
                       'LWNA0',
                       Node),
         Node = ifelse(Node == 'ICL',
                       'ICLB0',
                       Node),
         Node = ifelse(SiteID == 'NASONC',
                       'NALA0',
                       Node)) %>%
  distinct()


# Node network for DABOM
# site_df = writeTUMNodeNetwork()
site_df = writeTUMNodeNetwork_noUWE()


# Save file.
save(configuration,
     site_df,
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
yr = 2019
# start date is July 1 of the previous year
start_date = paste0(yr - 1, '0701')

# build parent-child table
par_ch_node = createParentChildDf(site_df,
                                  configuration,
                                  startDate = start_date)

root_node = par_ch_node %>%
  filter(nodeOrder == 1) %>%
  pull(ParentNode)

par_ch_site = par_ch_node %>%
  select(ParentNode, ChildNode) %>%
  mutate(ParentSite = str_remove(ParentNode, 'A0$'),
         ParentSite = str_remove(ParentSite, 'B0$'),
         ChildSite = str_remove(ChildNode, 'A0$'),
         ChildSite = str_remove(ChildSite, 'B0$')) %>%
  filter(ParentSite != ChildSite) %>%
  select(ParentSite,
         ChildSite) %>%
  distinct() %>%
  bind_rows(tibble(ParentSite = root_node,
                   ChildSite = root_node)) %>%
  left_join(configuration %>%
              select(ChildSite = SiteID,
                     SiteType,
                     RKM, RKMTotal,
                     lat = Latitude,
                     long = Longitude) %>%
              distinct()) %>%
  filter(!(ChildSite %in% ChildSite[duplicated(ChildSite)] & SiteType == 'MRR'))

# add nodes for black boxes
bbNodes = par_ch_site %>%
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
  bind_rows(bbNodes)

node_order = par_ch_site %>%
  rename(ParentNode = ParentSite,
         ChildNode = ChildSite) %>%
  getValidPaths(root_site = 'TUM') %>%
  createNodeOrder(configuration,
                  site_df,
                  step_num = 2) %>%
  select(-NodeSite) %>%
  mutate(Group = as.character(Group),
         Group = if_else(Node == 'TUM',
                         'TUM',
                         Group),
         Group = factor(Group,
                        levels = c('TUM', 'Peshastin', 'Icicle',
                                   'Chiwaukum', 'Chiwawa',
                                   'Nason', 'WhiteRiver', 'LittleWenatchee')),
         # Group = fct_explicit_na(Group),
         BranchNum = as.numeric(Group)) %>%
  arrange(Group, NodeOrder, RKM)

# build table of nodes
nodes = par_ch_site %>%
  select(ParentSite, ChildSite) %>%
  gather(type, node) %>%
  select(node) %>%
  distinct() %>%
  left_join(par_ch_site %>%
              select(-starts_with("Parent")) %>%
              rename(node = ChildSite)) %>%
  left_join(node_order %>%
              select(node = Node,
                     NodeOrder,
                     Group)) %>%
  mutate(NodeOrder = if_else(node == 'TUM_bb',
                         1.5,
                         as.numeric(NodeOrder))) %>%
  mutate(Group = if_else(node == 'TUM_bb',
                         'TUM',
                         as.character(Group))) %>%
  mutate(Group = factor(Group,
                        levels = levels(node_order$Group))) %>%
  arrange(Group, NodeOrder, RKM) %>%
  mutate(index = 1:n()) %>%
  rename(label = node) %>%
  select(index, label, everything())

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
# myLayout = c('partition')
# myLayout = c('circlepack')
myLayout = c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds',
             'randomly', 'fr', 'kk', 'drl', 'lgl')[4]

c('star', 'dh', 'lgl')

set.seed(8)
allGr_p = myGraph %>%
  # ggraph(layout = 'igraph',
  #        algorithm = 'nicely') +
  # ggraph(layout = 'igraph',
  #        algorithm = 'fr') +
  ggraph(layout = myLayout) +
  geom_edge_diagonal() +
  # geom_edge_link() +
  geom_node_point(aes(color = Group,
                      shape = nodeType),
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

#--------------------------------------------------
# netplot
set.seed(8)
# l = igraph::layout_with_fr(myGraph)
# l = igraph::layout_with_kk(myGraph)
# l = igraph::layout_with_dh(myGraph)
# l = igraph::layout_with_gem(myGraph)
# l = igraph::layout_nicely(myGraph)
# l = igraph::layout_on_grid(myGraph)
# l = igraph::layout_in_circle(myGraph)
# l = igraph::layout_with_graphopt(myGraph)
# l = igraph::layout_with_sugiyama(myGraph)
l = igraph::layout_as_tree(myGraph,
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

tum_sites = nplot(myGraph,
                  layout = l,
                  vertex.color = myColors[nodes$Group],
                  # vertex.color = myColors[nodes$nodeType],
                  # vertex.size = 1,
                  vertex.size.range = c(0.1, 0.05),
                  vertex.nsides = c(100, 3, 4, 4)[nodes$nodeType],
                  vertex.rot = c(0, 1.55, 0.78, 0.78)[nodes$nodeType],
                  vertex.label = nodes$label,

                             # vertex.label = myLabs,
                  vertex.label.fontsize = 12,
                  # edge.curvature = 0)
                  edge.curvature = pi/6)
tum_sites

pdf('analysis/figures/SpringChnk_TUM.pdf',
    width = 7,
    height = 6)
print(tum_sites)
dev.off()

# 0.Set environment ----------
#rm(list = ls())
#gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)


# 1.Figure D, F, G --------------------------------------------------------
## 1.1 load data -------
mycolors <- c(N = '#125919',
              L = '#446CDB',
              H = '#E0B131',
              C = '#CC2121',
              PC = '#E85ABD',
              CC = '#751E5B')

#protein matrix
dfpg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_demo_result/N20230320ProteomEx_v2_demo_whole_library_Yueliang/Whole_lib_report.pg_matrix.tsv')
colnames(dfpg) %<>% str_remove('^.+\\\\')

#protein info
protinfo <- dfpg %>% select(Protein.Group:First.Protein.Description)

#protein matrix (log2 transform)
mat <- dfpg %>%
  select(-(Protein.Ids:First.Protein.Description)) %>%
  column_to_rownames('Protein.Group') %>% log2()

#sample info
info <- rio::import('20230219_Batch_design_CRC_demo.xlsx', sheet = 'Sheet2') %>% select(-ncol(.))
info <- data.frame(SampleName = colnames(mat),
           Batch = str_extract(colnames(mat), 'b\\d+_(\\d+|pool)')) %>% 
  full_join(info)
info$BatchHead <- str_remove(info$Batch, '_.+$')
info$Slide %<>% factor()
info$Region %<>% factor(levels = c('N', 'L', 'H', 'C', 'PC', 'CC'), ordered = T)

#rename files in mat and regenerate dfpg
colnames(mat) %<>% str_extract('b\\d+_(\\d+|pool)')
dfpg <- mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .)


## 1.2 protein group numbers (D) ----
### boxplot ----
dfbar <- apply(mat, 2, function(y) sum(!is.na(y))) %>% as.data.frame() %>% 
  setNames('# proteins') %>% 
  rownames_to_column('Batch') %>% 
  inner_join(info, .) %>%
  filter(!is.na(Patient)) # remove pooled samples

#remove outliers in each region by their # proteins
get_outliers <- function(x){
  # x is a atomatic vector
  quant <- quantile(x)
  iqr <- IQR(x)
  c(lower = quant[2] - 1.5*iqr, upper = quant[4] + 1.5*iqr)
}
dfbar %<>% group_by(Region) %>% 
  summarise(lower = get_outliers(`# proteins`)[1],
            higher = get_outliers(`# proteins`)[2]) %>% 
  left_join(dfbar, .) %>% 
  mutate(Included = between(`# proteins`, lower, higher))

dfpg %<>% filter(Batch %in% dfbar$Batch[dfbar$Included])
mat <- mat[, dfpg$Batch]

lbls <- sort(unique(dfbar$Region))
dfcomp <- data.frame(X = lbls) %>% cross_join(data.frame(X = lbls)) %>%
  filter(X.x < X.y)
dfcomp$p.value <- dfcomp %>% 
  t() %>% as.data.frame() %>% as.list() %>% 
  sapply(function(x){
    x1 <- dfbar %>% filter(Region == x[1]) %>% pull(`# proteins`)
    x2 <- dfbar %>% filter(Region == x[2]) %>% pull(`# proteins`)
    t.test(x1, x2)$p.value
  })
comparisons <- dfcomp %>% filter(p.value < 0.05) %>% select(-p.value) %>% t() %>% as.data.frame() %>% as.list()

set.seed(10)
p4D <- ggplot(dfbar, aes(x = Region, y = `# proteins`, color = Region)) +
  stat_boxplot(geom = 'errorbar', width = 0.4)+
  geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
  geom_jitter(alpha = 0.7, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '# proteins') +
  scale_y_continuous(limits = c(0, 8200)) +
  scale_color_manual(values = mycolors) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('F4D_protein_number.pdf', p4D, width = 5, height = 5)  

tbl4 <- dfbar %>% cbind(NA) %>% cbind(NA) %>% 
  cbind(rbind(dfcomp, matrix(rep(NA, (nrow(dfbar) - nrow(dfcomp)) * ncol(dfcomp)), ncol = ncol(dfcomp)) %>% set_colnames(colnames(dfcomp))))
names(tbl4)[names(tbl4) == 'NA'] <- ''

### try heatmap style ----
library(pheatmap)

dfheat <- dfbar %>%
  filter(Included) %>% 
  select(Batch, Region, `# proteins`, Slide, Patient, BatchHead) %>%
  arrange(Region) %>%
  left_join(mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch')) %>% 
  set_rownames(.$Batch)
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

anno_colors <- list(
  Region = mycolors,
  BatchHead = brewer.pal(length(unique(dfheat$BatchHead)), 'Set2') %>% setNames(sort(unique(dfheat$BatchHead))),
  Patient = brewer.pal(length(unique(dfheat$Patient)), 'Set1') %>% setNames(sort(unique(dfheat$Patient))),
  '# proteins' = rev(brewer.pal(11, 'RdBu')),
  Slide = brewer.pal(nlevels(dfheat$Slide), 'Oranges') %>% setNames(levels(dfheat$Slide))
)

clust_col <- hclust(dist(t(mat_heat), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()

p4D_heat <- pheatmap(mat_heat,
         annotation_col = dfheat %>% select(`# proteins`:BatchHead, Region),
         annotation_colors = anno_colors,
         cluster_rows = T, cluster_cols = clust_col,
         clustering_distance_rows = 'euclidean',
         clustering_distance_cols = 'euclidean',
         clustering_method = 'ward.D2',
         show_rownames = F, show_colnames = T,
         color = c(brewer.pal(11,"RdYlBu")[9:7],brewer.pal(11,"RdYlBu")[4:2]),
         fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
         border_color = F,
         filename = 'F4D_whole_heatmap.pdf',
         width = 10, height = 10
)
tbl4_heat <- mat_heat[p4D_heat$tree_row$labels[p4D_heat$tree_row$order], p4D_heat$tree_col$labels[p4D_heat$tree_col$order]] %>% as.data.frame() %>% rownames_to_column()


## 1.3 DEA (F) ----
### missing ratio control -----
# single protein level control: <20% in at least 3 Regions
dfpg %>% count(Region)

dfpg_miss <- dfpg %>% group_by(Region) %>% select(-(SampleName:BatchHead)) %>%
  summarise_all(function(y) sum(is.na(y)) / length(y)) %>% 
  mutate(Region = as.character(Region))
dfpg_miss %<>% rbind(c(Region = 'ALL', apply(mat, 1, function(x) sum(is.na(x)) / length(x))))

tmp <- dfpg_miss %>%
  filter(Region != 'ALL') %>% 
  column_to_rownames('Region') %>%
  apply(2, function(y) sum(y < 0.2) >= 3 )
selected_pro <- names(tmp[tmp])
dfpg_miss %<>% select(Region, all_of(selected_pro))
mat <- mat[colnames(dfpg_miss)[-1], ]



# whole matrix level control
mr <- apply(mat, 1, function(x) sum(is.na(x)) / length(x))
plot(density(mr))

x <- seq(0.01, 0.99, 0.01)
y <- sapply(x, function(mr_cutoff){
  X <- mat[mr < mr_cutoff, ]
  return(sum(is.na(X)) / nrow(X) / ncol(X))
})
data.frame(x, y)
plot(x * 100, y * 100, xlab = '% NA ratio cutoff on protein level', ylab = '% whole matrix NA ratio')
abline(h = 6.5, v = 50)
# whole matrix NA ratio is 18.9% while the cutoff on protein level is set to 70%

mat50 <- mat[mr < 0.5, ]
dim(mat) # 4105  123
dim(mat50) # 4105  123

# #NA imputation with 0.8 * minimum * norm(1, 0.01)
# mat_naimp <- mat
# # set.seed(1)
# # mat_naimp[is.na(mat_naimp)] <- min(mat_naimp, na.rm = T) + log2(0.8) + log2(rnorm(sum(is.na(mat_naimp)), 1, 0.01))
# mat_naimp <- apply(mat, 1, function(x){
#   ifelse(is.na(x), min(x, na.rm = T) + log2(0.8), x)
# }) %>% t()
# dfpg_dea <- mat_naimp %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>%
#   inner_join(info, .) %>%
#   filter(!is.na(Region))

# NA ratio was considered by region respectively, and NA will be omitted
dfmiss <- dfpg %>%
  filter(!is.na(Region)) %>% 
  group_by(Region) %>% 
  select(-(SampleName:BatchHead)) %>% 
  summarise_all(function(y) sum(is.na(y)) / length(y))

dfpg_dea <- mat %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>%
  inner_join(info, .) %>%
  filter(!is.na(Region))


### DEA ------
# L/N, H/N, C/N
# L/N, H/L, C/H, CC/PC
comp_ls <- list(c('CC', 'N'), c('PC', 'N'), c('C', 'N'), c('H', 'N'), c('L', 'N'),
                c('H', 'L'), c('C', 'H'), c('PC', 'H'), c('CC', 'H'), c('CC', 'PC'))

DEA_ls <- lapply(comp_ls, function(comp){
  X <- dfpg_dea %>% filter(Region %in% comp) %>%
    select(-(SampleName:Slide), -(Rep:BatchHead)) %>% 
    pivot_longer(cols = -Region, names_to = 'Protein', values_to = 'Log2Intensity')
  ret <- plyr::ddply(X, 'Protein', function(dfsub){
    x1 <- dfsub$Log2Intensity[dfsub$Region == comp[1]]
    x2 <- dfsub$Log2Intensity[dfsub$Region == comp[2]]
    Log2FC <- log2(mean(2^x1, na.rm = T) / mean(2^x2, na.rm = T))
    P <- tryCatch(t.test(x1, x2)$p.value, error = function(e) NA)
    ret <- data.frame(Comparison = str_c(comp, collapse = ' vs. '),
                      Log2FC = Log2FC,
                      P = P)
    ret$adj.P <- p.adjust(ret$P, 'BH')
    return(ret)
  })
})
df_dea <- plyr::ldply(DEA_ls)
DEA <- protinfo %>%
  select(Protein.Ids, Genes) %>%
  rename(Protein = Protein.Ids) %>%
  inner_join(df_dea) %>%
  inner_join(dfmiss %>% column_to_rownames('Region') %>% t() %>% as.data.frame() %>% setNames(str_c('MissingRatio_', names(.))) %>% rownames_to_column('Protein'))
DEA %<>% mutate(
  Type = ifelse(Log2FC > 0, 'Up', 'Down'),
  is.significant = abs(Log2FC) > log2(2) & adj.P < 0.05,
  Label = str_c(Genes, ' (', Protein, ')'),
  Comparison = factor(Comparison, levels = sapply(comp_ls, str_c, collapse = ' vs. '), ordered = T)
)
DEA %>% filter(is.significant) %>% pull(Protein) %>% unique() %>% length() # 555 DEPs in total

DEP <- DEA %>% filter(is.significant)
DEP_solid <- apply(DEA, 1, function(x){
  if(is.na(x['Log2FC'])){
    ret <- FALSE
  } else if(x['Log2FC'] > 0){
    ret <- x[str_c('MissingRatio_', str_remove(x['Comparison'], ' vs\\..+$'))] < 0.5
  } else{
    ret <- x[str_c('MissingRatio_', str_remove(x['Comparison'], '^.+ vs\\. '))] < 0.5
  }
  return(ret)
}) %>% DEA[., ] %>% filter(is.significant)


### DEA visualization -----
#p value profile
ggplot(DEA) +
  aes(x = Comparison, y = adj.P, fill = Comparison) +
  geom_violin(alpha = 0.7, width = 0.5) +
  geom_hline(yintercept = 0.05, color = 'red') +
  scale_y_log10() +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic()

### DEP number ------
tmp <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[1:5]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)))) %>%
  select(Label, Comparison, Log2FC, adj.P, Type, is.significant) %>%
  count(Comparison, Type)
tbl6_line1 <- dfline1 <- tmp %>% group_by(Comparison) %>% summarise(nn = sum(n)) %>% 
  inner_join(tmp, .) %>% 
  mutate(x = str_c(Comparison, '\n(', nn, ')')) %>% 
  arrange(Comparison) %>% 
  mutate(x = factor(x, unique(x)))
pF1 <- ggplot(dfline1, aes(x = x, y = n, group = Type, color = Type))+
  geom_line(linewidth = 1) +
  geom_point(size = 4, alpha = 1) +
  labs(x = '', y = '', subtitle = '# differentially expressed proteins')+
  scale_color_manual(values = c(Up = 'red4', Down = 'blue4')) +
  theme_classic() +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 0, hjust = 0.5))

tmp <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[c(5:6, 8:10)]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)))) %>%
  select(Label, Comparison, Log2FC, adj.P, Type, is.significant) %>%
  count(Comparison, Type)
tbl6_line2 <- dfline2 <- tmp %>% group_by(Comparison) %>% summarise(nn = sum(n)) %>% 
  inner_join(tmp, .) %>% 
  mutate(x = str_c(Comparison, '\n(', nn, ')')) %>% 
  arrange(Comparison) %>% 
  mutate(x = factor(x, rev(unique(x))))
pF2 <- ggplot(dfline2, aes(x = x, y = n, group = Type, color = Type))+
  geom_line(linewidth = 1) +
  geom_point(size = 4, alpha = 1) +
  labs(x = '', y = '', subtitle = '# differentially expressed proteins')+
  scale_color_manual(values = c(Up = 'red4', Down = 'blue4')) +
  theme_classic() +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 0, hjust = 0.5))

pF <- ggpubr::ggarrange(pF1, pF2, nrow = 1, common.legend = T)
ggsave('F4F_DEPs_number.pdf', pF, width = 6, height = 3.5)


tbl6_line <- dfline <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)) %>%
  select(Label, Comparison, Log2FC, adj.P, Type, is.significant) %>%
  count(Comparison, Type)

PF_onePlot <- ggplot(dfline %>% filter(Comparison != 'C vs. H'), aes(x = Comparison, y = n, group = Type, color = Type))+
  geom_line(linewidth = 1) +
  geom_point(size = 4, alpha = 1) +
  labs(x = '', y = '', subtitle = '# differentially expressed proteins')+
  scale_color_manual(values = c(Up = 'red3', Down = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('F4F_DEPs_number_onePlot.pdf', PF_onePlot, width = 4, height = 3.5)

### heatmap -------
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

heat1 <- pheatmap(mat_heat[unique(DEP$Protein), ],
                  annotation_col = dfheat %>% select(Region),
                  annotation_colors = anno_colors,
                  cluster_rows = T, cluster_cols = F,
                  clustering_distance_rows = 'euclidean',
                  clustering_distance_cols = 'euclidean',
                  clustering_method = 'ward.D2',
                  show_rownames = F, show_colnames = F
)
heat2 <- pheatmap(mat_heat[unique(DEP_solid$Protein), ],
                  annotation_col = dfheat %>% select(Region),
                  annotation_colors = anno_colors,
                  cluster_rows = T, cluster_cols = F,
                  clustering_distance_rows = 'euclidean',
                  clustering_distance_cols = 'euclidean',
                  clustering_method = 'ward.D2',
                  show_rownames = F, show_colnames = F
)
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
# mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlGn")))(50)
bk <- unique(c(seq(-2, 2, length = 50)))
pF_heat1 <- pheatmap(mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ],
                     annotation_col = dfheat %>% select(Region),
                     annotation_colors = anno_colors,
                     cluster_rows = F, cluster_cols = F,
                     show_rownames = F, show_colnames = F,
                     color = heat_colors, breaks = bk,
                     fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                     border_color = F, na_col = '#AAAAAA',
                     gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                     filename = 'F4F_555DEPs_heatmap.pdf',
                     width = 10, height = 5
)


clust_col <- hclust(dist(t(mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ]), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
pF_heat1_cls <- pheatmap(
  mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ],
  scale = 'none',
  annotation_col = dfheat %>% select(Region, Patient),
  annotation_colors = anno_colors,
  cluster_rows = T, cluster_cols = clust_col,
  clustering_method = 'ward.D2',
  clustering_distance_cols = 'euclidean', clustering_distance_rows = 'euclidean',
  show_rownames = F, show_colnames = F,
  # color = heat_colors, breaks = bk,
  fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
  border_color = F, na_col = '#AAAAAA',
  filename = 'F4F_555DEPs_heatmap_hclust.pdf',
  width = 10, height = 5)

  
p4F_heat2 <- pheatmap(mat_heat[heat2$tree_row$labels[heat2$tree_row$order], ],
                     annotation_col = dfheat %>% select(Region),
                     annotation_colors = anno_colors,
                     cluster_rows = F, cluster_cols = F,
                     show_rownames = F, show_colnames = F,
                     color = heat_colors, breaks = bk,
                     fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                     border_color = F, na_col = '#AAAAAA',
                     gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                     filename = 'F4F_516of555DEPs_heatmap.pdf',
                     width = 10, height = 5
)

#### UMAP ------
library(umap)
X_imp <- mat_heat[heat1$tree_row$labels[heat1$tree_row$order], ]
X_imp[is.na(X_imp)] <- 0

umap <- umap(t(X_imp), n_neighbors = 56, min_dist = 0.5) # default n_neighbors = 10, min_dist = 0.1
df_umap <- info %>% inner_join(umap$layout %>% as.data.frame() %>% rownames_to_column('Batch'))
df_umap_center <- df_umap %>% group_by(Region) %>% summarise_at(vars(center_x = V1, center_y = V2), mean) # calculate center of clusters
df_umap %<>% left_join(df_umap_center)

pF_umap <- ggplot(df_umap, aes(x = V1, y = V2, color = Region)) +
  geom_point(size = 4, alpha = 0.9)+
  # geom_point(data = df_umap_center, aes(x = center_x, y = center_y),
  #            shape = 'square', size = 4, alpha = 0.5, show.legend = F) +
  # geom_curve(data = df_umap,
  #            aes(xend = V1, yend = V2, x = center_x, y = center_y),
  #            alpha = 0.6, show.legend = F, curvature = 0.5)+
  # stat_ellipse(type = 'norm', level = 0.95, size = 1)+
  labs(x = 'UMAP1', y = 'UMAP2') +
  scale_color_manual(values = mycolors) +
  theme_classic() +
  theme(text = element_text(size = 10), legend.position = 'none')
ggsave('F4F_555DEPs_UMAP.pdf', pF_umap, width = 4, height = 4)

pF_umap_pat <- ggplot(df_umap, aes(x = V1, y = V2, color = Patient)) +
  geom_point(size = 4, alpha = 0.9)+
  # geom_point(data = df_umap_center, aes(x = center_x, y = center_y),
  #            shape = 'square', size = 4, alpha = 0.5, show.legend = F) +
  # geom_curve(data = df_umap,
  #            aes(xend = V1, yend = V2, x = center_x, y = center_y),
  #            alpha = 0.6, show.legend = F, curvature = 0.5)+
  # stat_ellipse(type = 'norm', level = 0.95, size = 1)+
  labs(x = 'UMAP1', y = 'UMAP2') +
  # scale_color_manual(values = mycolors) +
  theme_classic() +
  theme(text = element_text(size = 10), legend.position = 'right')
ggsave('F4F_555DEPs_UMAP_patient.pdf', pF_umap_pat, width = 4.2, height = 4)



### other figures -----
DEA_wide <- DEA %>% select(Label, Comparison, Log2FC, adj.P) %>% 
  pivot_wider(id_cols = Label, names_from = Comparison, values_from = c('Log2FC', 'adj.P'))
DEA_wide %>% 
  select(Label, matches('^adj\\.P')) %>% filter_if(is.double, function(x) x < 0.05) %>% 
  semi_join(DEA_wide, .) %>% 
  select(Label, matches('^Log2FC')) %>% filter_if(is.double, function(x) abs(x) > log2(1.5))



set.seed(10)
pF_volcano1 <- DEA %>% 
  filter(Comparison %in% levels(DEA$Comparison)[c(1:5)]) %>% 
  mutate(Comparison = factor(Comparison, rev(levels(DEA$Comparison)))) %>% 
  ggplot(aes(x = Comparison, y = Log2FC, color = interaction(Type, is.significant)))+
  # geom_jitter(alpha = 0.5, size = 0.3, width = 0.4) +
  geom_jitter(alpha = 0.6, size = 3, width = 0.4) +
  labs(x = '', y = 'Log2 Fold-Change', subtitle = '# differentially expressed proteins')+
  scale_color_manual(name = '',
                     values = c(Up.FALSE = '#AAAAAA', Down.FALSE = '#AAAAAA',
                                Up.TRUE = 'red4', Down.TRUE = 'blue4'),
                     labels = c('', '|FC|<2 or adj.P>0.05', 'Down', 'Up')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'top') #+
# geom_jitter(
#   data = DEP %>% filter(
#     Comparison %in% levels(DEA$Comparison)[c(1:5)],
#     Protein %in% c('P06731', 'Q01628', 'P25815', 'A8K7I4', 'Q9Y6R7')
#   ),
#   size = 2, width = 0.4) +
# ggrepel::geom_text_repel(
#   data = DEP %>% filter(
#     Comparison %in% levels(DEA$Comparison)[c(1:5)],
#     Protein %in% c('P06731', 'Q01628', 'P25815', 'A8K7I4', 'Q9Y6R7')
#   ),
#   aes(label = Label), color = 'black', size = 3,
#   force = 1.2, seed = 10,
#   arrow = arrow(length = unit(0.008, "npc"),
#                 type = "open", ends = "last")
# )
set.seed(10)
pF_volcano2 <- DEA %>% 
  filter(Comparison %in% levels(DEA$Comparison)[c(5:6, 8:10)]) %>% 
  ggplot(aes(x = Comparison, y = Log2FC, color = interaction(Type, is.significant)))+
  # geom_jitter(alpha = 0.5, size = 0.3, width = 0.4) +
  geom_jitter(alpha = 0.6, size = 3, width = 0.4) +
  labs(x = '', y = 'Log2 Fold-Change', subtitle = '# differentially expressed proteins')+
  scale_color_manual(name = '',
                     values = c(Up.FALSE = '#AAAAAA', Down.FALSE = '#AAAAAA',
                                Up.TRUE = 'red4', Down.TRUE = 'blue4'),
                     labels = c('', '|FC|<2 or adj.P>0.05', 'Down', 'Up')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'top') #+
# geom_jitter(
#   data = DEP %>% filter(
#     Comparison %in% levels(DEA$Comparison)[c(1:5)],
#     Protein %in% c('P06731', 'Q01628', 'P25815', 'A8K7I4', 'Q9Y6R7')
#   ),
#   size = 2, width = 0.4) +
# ggrepel::geom_text_repel(
#   data = DEP %>% filter(
#     Comparison %in% levels(DEA$Comparison)[c(1:5)],
#     Protein %in% c('P06731', 'Q01628', 'P25815', 'A8K7I4', 'Q9Y6R7')
#   ),
#   aes(label = Label), color = 'black', size = 3,
#   force = 1.2, seed = 10,
#   arrow = arrow(length = unit(0.008, "npc"),
#                 type = "open", ends = "last")
# )
pF_volcano <- ggpubr::ggarrange(pF_volcano1, pF_volcano2, nrow = 1, ncol = 2, common.legend = T)
ggsave('F4F_DEPs_volcano.pdf', pF_volcano, width = 12, height = 6)

library(patchwork)
pF_volcano_line <- wrap_plots(pF_volcano, pF) + plot_layout(heights = c(3.5,1))
ggsave('F4F_DEPs_volcano_lines.pdf', pF_volcano_line, width = 12, height = 7.7)


# dfpg_dea$Region
# mat_naimp[, dfpg_dea$Batch[dfpg_dea$Region == comp[1]]]

## 1.4 upset of DEPs (G) ----
tmp <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[c(1:10)]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)[c(1:10)]))) %>% 
  select(Protein, Comparison) %>%
  distinct()
DEP_list <- split(tmp$Protein, tmp$Comparison)
names(DEP_list) <- str_c(names(DEP_list), ' (', sapply(DEP_list, length), ')')
pdf('F4G_DEPs_upset.pdf', width = 10, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(DEP_list),
    sets = names(DEP_list), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on',
    order.by = 'freq',# order.by = c('freq', 'degree'),
    mainbar.y.label = '# proteins',
    point.size = 3, 
    line.size = 1,
    set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5),
    mb.ratio = c(0.6, 0.4),
    queries = list(
      list(query = UpSetR::intersects, params = names(DEP_list)[6:10], color = "#F39800", font.color = "#F39800", active = T, query.name = 'a')
      # list(query = UpSetR::intersects, params = names(DEP_list)[3:6], color = "#F39800", font.color = "#F39800", active = T, query.name = 'b')
    ),
    # query.legend = "none",
  )
)
graphics.off()

tmp1 <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[c(1:5)]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)[c(1:5)]))) %>% 
  select(Protein, Comparison) %>%
  distinct()
DEP_list1 <- split(tmp1$Protein, tmp1$Comparison)
names(DEP_list1) <- str_c(names(DEP_list1), ' (', sapply(DEP_list1, length), ')')

tmp2 <- DEP %>%
  filter(Comparison %in% levels(DEA$Comparison)[c(5:10)]) %>%
  mutate(Comparison = factor(Comparison, levels = rev(levels(DEA$Comparison)[c(5:10)]))) %>% 
  select(Protein, Comparison) %>%
  distinct()
DEP_list2 <- split(tmp2$Protein, tmp2$Comparison)
names(DEP_list2) <- str_c(names(DEP_list2), ' (', sapply(DEP_list2, length), ')')

pdf('F4G_DEPs_upset_split.pdf', width = 10, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(DEP_list1),
    sets = rev(names(DEP_list1)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on',
    order.by = 'freq',# order.by = c('freq', 'degree'),
    mainbar.y.label = '# proteins',
    point.size = 3, 
    line.size = 1,
    set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5),
    mb.ratio = c(0.6, 0.4),
    # queries = list(
    #   list(query = UpSetR::intersects, params = names(DEP_list)[1:5], color = "#F39800", font.color = "#F39800", active = T, query.name = 'a')
    #   # list(query = UpSetR::intersects, params = names(DEP_list)[3:6], color = "#F39800", font.color = "#F39800", active = T, query.name = 'b')
    # ),
    # query.legend = "none",
  )
)
print(
  UpSetR::upset(
    UpSetR::fromList(DEP_list2),
    sets = rev(names(DEP_list2)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on',
    order.by = 'freq',# order.by = c('freq', 'degree'),
    mainbar.y.label = '# proteins',
    point.size = 3, 
    line.size = 1,
    set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5),
    mb.ratio = c(0.6, 0.4)
  )
)
graphics.off()




list(
  Figure4D = tbl4,
  Figure4D_heatmap = tbl4_heat,
  Figure4F_progress = tbl6_line1,
  Figure4F_othersVSn = tbl6_line2,
  Figure4F_DEA = DEA,
  Figure4F_DEP = DEP,
  Figure4F_heatmap = mat_heat[pF_heat1_cls$tree_row$labels[pF_heat1_cls$tree_row$order],
                              pF_heat1_cls$tree_col$labels[pF_heat1_cls$tree_col$order]] %>% 
    as.data.frame() %>% 
    rownames_to_column('Protein') %>% inner_join(DEP %>% distinct(Protein, Genes), .),
  Figure4G = tmp %>% mutate(value = T) %>% pivot_wider(names_from = 'Comparison', values_from = 'value'),
  Figure4G_progress = tmp1 %>% mutate(value = T) %>% pivot_wider(names_from = 'Comparison', values_from = 'value'),
  Figure4G_othersVSn = tmp2 %>% mutate(value = T) %>% pivot_wider(names_from = 'Comparison', values_from = 'value')
) %>%
  rio::export('F4_D_F_G.xlsx')

save.image('Figure4.RData')




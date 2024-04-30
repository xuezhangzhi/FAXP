# 0.Set environment ----------
#rm(list = ls())
#gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)




# 1.DEA -----------------------------------------------------
load('../Figure4/Figure4.RData') # based on Figure4 data

## 1.1 statistics ------------------
### figure4 results ----
dim(DEA) # 41050 15
DEP %>% pull(Protein) %>% unique() %>% length() # 555 DEPs in total
DEP_solid %>% pull(Protein) %>% unique() %>% length() # 516

### disease progress associated proteins -----
my_comparisons <- c('L vs. N', 'H vs. L', 'PC vs. H', 'CC vs. H')
DEA_prog <- DEA %>% filter(Comparison %in% my_comparisons,
                           adj.P < 0.05)
prot_prog <- DEA_prog %>%
  pivot_wider(id_cols = 'Protein', names_from = 'Comparison', values_from = 'Log2FC') %>% 
  column_to_rownames('Protein') %>% 
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>% # do not consider NA fold-change
  apply(1, function(x) (all(x >= 0) | all(x <= 0)) & any(abs(x) > log2(2))) %>% # fold-change values should be in the same trend, and at least one of them >2
  which() %>% names()
DEP_prog <- DEA_prog %>% filter(Protein %in% prot_prog)


### boxplot ------
prot_prog1 <- DEP_prog %>% count(Protein) %>% filter(n > 1) %>% pull(Protein) # at least two significant comparison
DEP_prog1 <- DEP_prog %>% filter(Protein %in% prot_prog1)

prot_prog2 <- DEP_prog1 %>% select(matches('^Missing.+_')) %>% apply(1, function(x) sum(x < 0.2) >= 3) %>% filter(DEP_prog1, .) %>% pull(Protein) %>% unique() # Missing <20% in at least three Regions
DEP_prog2 <- DEP_prog1 %>% filter(Protein %in% prot_prog2)

dfbox <- dfpg %>%
  filter(
    !is.na(Region), Region != 'C'
  ) %>% 
  select(SampleName:BatchHead, all_of(prot_prog2)) %>% 
  pivot_longer(cols = -c(SampleName:BatchHead), names_to = 'Protein', values_to = 'Log2Intensity') %>% 
  left_join(DEA %>% distinct(Protein, Genes, Label), relationship = 'many-to-many') %>%
  mutate(Protein = factor(Protein, levels = prot_prog2))

interested_genes <- c('CEACAM5','CEACAM6','FCGBP','H1-1','S100P','S100A4')
nafill <- dfpg %>% select(-(SampleName:BatchHead)) %>%
  min(na.rm = T) %>% `+`(log2(0.8))
boxplot_ls_main <- dfbox %>% 
  filter(Genes %in% interested_genes) %>% 
  mutate(Genes = factor(Genes, interested_genes)) %>% 
  arrange(Genes) %>% 
  mutate(Protein = factor(Protein, unique(Protein))) %>% 
  plyr::dlply('Protein', function(dfsub){
    dfna <- dfsub %>% filter(is.na(Log2Intensity)) %>% 
      mutate(Log2Intensity = nafill)
    ymax <- max(dfsub$Log2Intensity, na.rm = T)
    DEPsub <- DEA %>%
      filter(
        Protein == dfsub$Protein[1], adj.P < 0.05,
        Comparison %in% c('L vs. N', 'H vs. L', 'PC vs. H', 'CC vs. H',
                          'CC vs. PC')
      ) %>%
      mutate(y = ymax * seq(1.05, by = 0.08, length.out = nrow(.)))
    label <- dfsub$Label[1]
    set.seed(10)
    ggplot(dfsub, aes(x = Region, y = `Log2Intensity`, color = Region)) +
      stat_boxplot(geom = 'errorbar', width = 0.4)+
      geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
      geom_jitter(alpha = 0.4, size = 1, width = 0.15)+
      geom_jitter(data = dfna, alpha = 0.4, size = 1, width = 0.15) +
      labs(x = '', y = '', subtitle = label) +
      scale_y_continuous(limits = c(0, ymax * 1.3)) +
      scale_color_manual(values = mycolors) +
      theme_classic() +
      theme(text = element_text(size = 15), legend.position = 'none') +
      # ggpubr::stat_compare_means(
      #   method = 't.test',
      #   map_signif_level = F, hide.ns = F, label = 'p.signif',
      #   comparisons = str_split(my_comparisons, ' vs. '),
      #   size = 5, hjust = 0.5, vjust = 0) +
      ggsignif::geom_signif(
        data = DEPsub,
        manual = T, inherit.aes = F,
        aes(xmax = str_remove(Comparison, '^.+vs\\. '),
            xmin = str_remove(Comparison, ' vs\\. .+$'),
            y_position = y,
            annotations = sprintf('%.2e', adj.P))
      )
  })
boxplots <- ggpubr::ggarrange(plotlist = boxplot_ls_main, nrow = 2, ncol = 3)
ggsave('F5A_DEPs_progress_box_main.pdf', boxplots, width = 3 * 4, height = 2 * 4, limitsize = F)

boxplot_ls_supp <- dfbox %>% 
  filter(!(Genes %in% interested_genes)) %>% 
  mutate(Genes = factor(Genes, interested_genes)) %>% 
  arrange(Genes) %>% 
  mutate(Protein = factor(Protein, unique(Protein))) %>% 
  plyr::dlply('Protein', function(dfsub){
    dfna <- dfsub %>% filter(is.na(Log2Intensity)) %>% 
      mutate(Log2Intensity = nafill)
    ymax <- max(dfsub$Log2Intensity, na.rm = T)
    DEPsub <- DEA %>%
      filter(
        Protein == dfsub$Protein[1], adj.P < 0.05,
        Comparison %in% c('L vs. N', 'H vs. L', 'PC vs. H', 'CC vs. H',
                          'CC vs. PC')
      ) %>%
      mutate(y = ymax * seq(1.05, by = 0.08, length.out = nrow(.)))
    label <- dfsub$Label[1]
    set.seed(10)
    ggplot(dfsub, aes(x = Region, y = `Log2Intensity`, color = Region)) +
      stat_boxplot(geom = 'errorbar', width = 0.4)+
      geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
      geom_jitter(alpha = 0.4, size = 1, width = 0.15)+
      geom_jitter(data = dfna, alpha = 0.4, size = 1, width = 0.15) +
      labs(x = '', y = '', subtitle = label) +
      scale_y_continuous(limits = c(0, ymax * 1.3)) +
      scale_color_manual(values = mycolors) +
      theme_classic() +
      theme(text = element_text(size = 15), legend.position = 'none') +
      # ggpubr::stat_compare_means(
      #   method = 't.test',
      #   map_signif_level = F, hide.ns = F, label = 'p.signif',
      #   comparisons = str_split(my_comparisons, ' vs. '),
      #   size = 5, hjust = 0.5, vjust = 0) +
      ggsignif::geom_signif(
        data = DEPsub,
        manual = T, inherit.aes = F,
        aes(xmax = str_remove(Comparison, '^.+vs\\. '),
            xmin = str_remove(Comparison, ' vs\\. .+$'),
            y_position = y,
            annotations = sprintf('%.2e', adj.P))
      )
  })
boxplots_supp <- ggpubr::ggarrange(plotlist = boxplot_ls_supp, nrow = 9, ncol = 5)
ggsave('F5A_DEPs_progress_box_supp.pdf', boxplots_supp, width = 5 * 4, height = 9 * 4, limitsize = F)


# ### add median value filtering
# dfbox %>%
#   group_by(Protein, Region) %>% 
#   summarise_at(vars(Log2MedianIntensity = Log2Intensity), median, na.rm = T) %>% 
#   group_by(Protein) %>%
#   summarise_at(vars(Mean_Median_both = Log2MedianIntensity), function(x) identical(x, sort(x))) %>% 
#   filter(Mean_Median_both)
# # Q9BQL6 only
# 
# 
# interested_labels <- c('P06731', 'Q01628', 'P25815', 'A8K7I4', 'Q9Y6R7')
# intersect(interested_labels, prot_prog) # not include A8K7I4
# boxplot_ls <- plyr::dlply(dfbox %>% filter(Protein %in% intersect(interested_labels, prot_prog)), 'Protein', function(dfsub){
#   label <- dfsub$Label[1]
#   set.seed(10)
#   ggplot(dfsub, aes(x = Region, y = `Log2Intensity`, color = Region)) +
#     stat_boxplot(geom = 'errorbar', width = 0.4)+
#     geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
#     geom_jitter(alpha = 0.4, size = 1, width = 0.15)+
#     labs(x = '', y = '', subtitle = label) +
#     scale_color_manual(values = mycolors) +
#     theme_classic() +
#     theme(text = element_text(size = 15), legend.position = 'none') +
#     ggpubr::stat_compare_means(
#       method = 't.test',
#       map_signif_level = F,# hide.ns = F, label = 'p.signif',
#       comparisons = str_split(my_comparisons, ' vs. '),
#       size = 5, hjust = 0.5, vjust = 0)
# })
# boxplots <- ggpubr::ggarrange(plotlist = boxplot_ls, nrow = 2, ncol = 2)
# ggsave('F5A_DEPs_progress_box_CEACAM5_S100P_IFITM3_FCGBP.pdf', boxplots, width = 2 * 4, height = 2 * 4)

## 1.2 volcano ------------
# # DEA_wide <- DEP_prog2 %>% select(Label, Comparison, Log2FC, adj.P) %>% 
# #   pivot_wider(id_cols = Label, names_from = Comparison, values_from = c('Log2FC', 'adj.P'))
# # DEA_wide %>% 
# #   select(Label, matches('^adj\\.P')) %>% filter_if(is.double, function(x) x < 0.05) %>% 
# #   semi_join(DEA_wide, .) %>% 
# #   select(Label, matches('^Log2FC')) %>% filter_if(is.double, function(x) abs(x) > log2(1.2))
# 
# set.seed(10)
# p5A_volcano <- DEA %>% filter(Comparison %in% my_comparisons) %>% 
#   ggplot(aes(x = Comparison, y = Log2FC, color = interaction(Type, is.significant)))+
#   geom_jitter(#data = DEA %>% filter(abs(Log2FC) > 2),
#     alpha = 0.5, size = 0.5, width = 0.4) +
#   labs(x = '', y = 'Log2 Fold-Change', subtitle = '# differentially expressed proteins')+
#   scale_color_manual(name = '',
#                      values = c(Up.FALSE = '#DDDDDD', Down.FALSE = '#DDDDDD',
#                                 Up.TRUE = 'red3', Down.TRUE = 'blue3'),
#                      labels = c('', '|FC|<2 or adj.P>0.05', 'Down', 'Up')) +
#   theme_classic() +
#   theme(text = element_text(size = 15), legend.position = 'top') +
#   # geom_jitter(data = DEP_prog2 %>% filter(is.significant),
#   #             size = 3, width = 0.4) +
#   ggrepel::geom_text_repel(
#     data = DEP_prog2 %>% filter(is.significant),
#     aes(label = Genes, size = abs(Log2FC)),
#     direction = 'x', show.legend = F,
#     force = 1.2, seed = 10, max.overlaps = 100, segment.color = NA
#   ) +
#   scale_size_continuous(range = c(0.8, 4), breaks = c(0.5, 3))
# ggsave('F5A_DEPs_progress_volcano.pdf', p5A_volcano, width = 10, height = 8)


## 1.3 heatmap ---------------------
library(pheatmap)
mat_heat <- dfheat %>%
  filter(Region != 'C') %>% 
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

clust_col <- hclust(dist(t(mat_heat[unique(DEP_prog2$Protein), ]), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat %>% filter(Region != 'C'))) %>% # reorder to keep original order as much as possible
  as.hclust()
heat3 <- pheatmap(mat_heat[unique(DEP_prog2$Protein), ],
                  annotation_col = dfheat %>% filter(Region != 'C') %>% select(Region),
                  annotation_colors = anno_colors,
                  cluster_rows = T, cluster_cols = clust_col,
                  clustering_distance_rows = 'euclidean',
                  clustering_distance_cols = 'euclidean',
                  clustering_method = 'ward.D2',
                  show_rownames = F, show_colnames = F
)
mat_heat <- dfheat %>%
  filter(Region != 'C') %>% 
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
# mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

heat_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdYlGn")))(50)
bk <- unique(c(seq(-2, 2, length = 50)))
X <- mat_heat[heat3$tree_row$labels[heat3$tree_row$order], ]
rownames(X) <- DEP_prog2 %>% distinct(Protein, Label) %>% set_rownames(.$Protein) %>% .[rownames(X), ] %>% pull(Label)
pF_heat3 <- pheatmap(X,
                     annotation_col = dfheat %>% filter(Region != 'C') %>% select(Region, Patient),
                     annotation_colors = anno_colors,
                     cluster_rows = F, cluster_cols = F,
                     show_rownames = T, show_colnames = F,
                     color = heat_colors, breaks = bk,
                     fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                     border_color = F, na_col = '#AAAAAA',
                     gaps_col = cumsum(table(dfheat$Region[dfheat$Region != 'C']))[-length(cumsum(table(dfheat$Region)))],
                     filename = 'F5A_DEPs_progress_heatmap.pdf',
                     width = 10, height = 8
)

pF_heat3 <- pheatmap(X,
                     annotation_col = dfheat %>% filter(Region != 'C') %>% select(Region, Patient),
                     annotation_colors = anno_colors,
                     cluster_rows = F, cluster_cols = clust_col,
                     show_rownames = T, show_colnames = F,
                     color = heat_colors, breaks = bk,
                     fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                     border_color = F, na_col = '#AAAAAA',
                     filename = 'F5A_DEPs_progress_hclust_heatmap.pdf',
                     width = 10, height = 8
)

## 1.4 UMAP ---------
library(umap)
X_imp <- X
X_imp[is.na(X_imp)] <- 0

umap <- umap(t(X_imp), n_neighbors = 50, min_dist = 0.1) # default n_neighbors = 10, min_dist = 0.1
df_umap <- info %>% inner_join(umap$layout %>% as.data.frame() %>% rownames_to_column('Batch'))
df_umap_center <- df_umap %>% group_by(Region) %>% summarise_at(vars(center_x = V1, center_y = V2), mean) # calculate center of clusters
df_umap %<>% left_join(df_umap_center)

p_umap <- ggplot(df_umap, aes(x = V1, y = V2, color = Region)) +
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
ggsave('F5A_DEPs_progress_UMAP.pdf', p_umap, width = 4, height = 4)

p_umap <- ggplot(df_umap, aes(x = V1, y = V2, color = Patient)) +
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
  theme(text = element_text(size = 10))
ggsave('F5A_DEPs_progress_UMAP_patient.pdf', p_umap, width = 4.2, height = 4)

## 1.5 paper query ---------
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(UniprotR)
library(easyPubMed)

# get Protein Functions from UNIPROT database
uniprot_pf <- GetProteinFunction(prot_prog2)
uniprot_pf$Function..CC.

# Publication
uniprot_pub <- GetPublication(prot_prog2, directorypath = './')

reticulate::source_python('data/article_from_pmid.py') # use python function
pubmed_articles <- list()
for(i in 1:nrow(uniprot_pub)){
  cat(i, '...\r')
  pmids <- str_split(uniprot_pub$ProteinDataTable[i], '; ')[[1]]
  pubmed_articles[[i]] <- article_from_pmid(pmids)
}
df_articles <- pubmed_articles %>% setNames(rownames(uniprot_pub)) %>% 
  plyr::ldply(.id = 'Protein')

df_articles[, 'Keyword: Tumor Progress'] <- str_detect(df_articles$Abstracts, '([Cc]ancer|[Tt]umor) progress')
df_articles[, 'Keyword: Tumor Development'] <- str_detect(df_articles$Abstracts, '([Cc]ancer|[Tt]umor) development')
df_articles[, 'Keyword: Tumor Invasion'] <- str_detect(df_articles$Abstracts, '([Cc]ancer|[Tt]umor) invasion')
df_articles[, 'Keyword: Tumor Metastasis'] <- str_detect(df_articles$Abstracts, '([Cc]ancer|[Tt]umor) metastasis')
df_articles[, 'Keyword: Cell growth'] <- str_detect(df_articles$Abstracts, '[Cc]ell growth')
df_articles[, 'Keyword: Cell Adhesion'] <- str_detect(df_articles$Abstracts, '[Cc]ell adhesion')
df_articles[, 'Keyword: Cell Migration'] <- str_detect(df_articles$Abstracts, '[Cc]ell migration')

df_articles <- DEP_prog2 %>% filter(adj.P < 0.05) %>%
  select(Protein, Genes, Comparison, Log2FC) %>%
  rename(Gene = Genes) %>% 
  pivot_wider(id_cols = c('Protein', 'Gene'), names_from = 'Comparison', values_from = 'Log2FC') %>% 
  left_join(df_articles, .)

df_articles1 <- df_articles %>%
  select(starts_with('Keyword')) %>%
  apply(1, function(x) any(x)) %>% 
  df_articles[., ]

length(unique(df_articles$PMID))
length(unique(df_articles1$PMID))
list(article744 = df_articles,
     article34_progressKeyword = df_articles1,
     UNIPROT_annotation = uniprot_pf %>% select(Function..CC.) %>% rownames_to_column('Protein') %>% left_join(df_articles %>% distinct(Protein, Gene), .)) %>% 
  rio::export('F5A_DEP_disease_progress_proteins50_getArticles.xlsx')



## 1.6 Pathway enrichment -----------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

### Metascape -------
writeClipboard(unique(DEP_prog2$Genes))
writeClipboard(unique(DEP_prog$Genes))

### GO ------------
# Gene ontology profile
GO <- enrichGO(unique(DEP_prog$Genes), OrgDb = org.Hs.eg.db,
               ont = 'ALL', keyType = 'SYMBOL', readable = T,
               minGSSize = 5, maxGSSize = 5000, pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GO_simplify <- clusterProfiler::simplify(GO, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgo <- GO[]
dfgo_simplify <- GO_simplify[]
enrichplot::dotplot(GO, showCategory = 5, split = 'ONTOLOGY', title = 'Gene Ontology') + facet_grid(ONTOLOGY~., scale='free')

# focus on GOBP
GOBP <- enrichGO(unique(DEP_prog2$Genes), OrgDb = org.Hs.eg.db,
                 ont = 'BP', keyType = 'SYMBOL', readable = T,
                 minGSSize = 10, maxGSSize = 2500, pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOBP_simplify <- clusterProfiler::simplify(GOBP, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgobp <- GOBP[]
dfgobp_simplify <- GOBP_simplify[]

plot_gobp <- enrichplot::dotplot(GOBP_simplify, showCategory = 8, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'GOBP (simplified)') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
ggsave('F5A_DEPs_progress_109_GOBP.pdf', plot_gobp, width = 4.5, height = 3.5)

pdf('F5A_DEPs_progress_109_GOBP_similar.pdf', width = 7, height = 4)
similar_matrix <- simplifyEnrichment::GO_similarity(dfgobp_simplify$ID, ont = 'BP', db = 'org.Hs.eg.db')
GOBP_simplify_sim <- simplifyEnrichment::simplifyGO(similar_matrix)
graphics.off()

### KEGG ------------
library(KEGG.db)

tmp <- clusterProfiler::bitr(DEP_prog$Genes, OrgDb = org.Hs.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID')
# tmp <- clusterProfiler::bitr_kegg(DEP_prog$Protein, fromType = 'uniprot', toType = 'ncbi-proteinid', organism = 'hsa')

# getOption("clusterProfiler.download.method")
# R.utils::setOption("clusterProfiler.download.method","wininet")
# R.utils::setOption("clusterProfiler.download.method","auto")
# R.utils::setOption("clusterProfiler.download.method","wget")
# R.utils::setOption("clusterProfiler.download.method","libcurl")
# R.utils::setOption("clusterProfiler.download.method","curl")
# R.utils::setOption("clusterProfiler.download.method","internal")
KEGG <- enrichKEGG(tmp$ENTREZID, organism = 'hsa', use_internal_data = T,
                   keyType = 'ncbi-geneid',
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dfkegg <- KEGG[]

plot_kegg <- enrichplot::dotplot(KEGG, showCategory = 5, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'KEGG') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
ggsave('F5A_DEPs_progress_109_KEGG.pdf', plot_kegg, width = 4.5, height = 3)



# # check hsa05210 (Colorectal cancer - Homo sapiens (human))
# library(limma)
# tab <- getGeneKEGGLinks(species = "hsa")
# tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,
#                      column = "SYMBOL", keytype = "ENTREZID")
# tab_ <- tab %>% full_join(DEA, by = c(Symbol = 'Genes'), relationship = 'many-to-many')
# tab_ %>% dplyr::filter(PathwayID == 'hsa05210') %>% View()


### Reactome ------
library(ReactomePA)
# gsePathway()

tmp <- clusterProfiler::bitr(DEP_prog2$Genes, OrgDb = org.Hs.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID')
Reactome <- enrichPathway(tmp$ENTREZID,
                          organism = 'human', readable = T,
                          minGSSize = 10, maxGSSize = 2500, pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dfreactome <- Reactome[]

plot_reactome <- enrichplot::dotplot(Reactome, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'Reactome') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
ggsave('F5A_DEPs_progress_50_Reactome.pdf', plot_reactome, width = 5, height = 5)

# view interested pathways
gList <- DEP_prog %>%
  filter(abs(Log2FC) > log2(2)) %>% 
  left_join(tmp, c(Genes = 'SYMBOL'), relationship = 'many-to-many') %>% 
  distinct(ENTREZID, Log2FC) %>%
  drop_na() %>%
  arrange(Log2FC) %>%
  distinct(ENTREZID, .keep_all = T)

viewPathway("Extracellular matrix organization",
            readable = TRUE,
            foldChange = gList$Log2FC %>% setNames(gList$ENTREZID))

pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)


## 1.7 GSVA ------
# prepare gene matrix
dfgg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_demo_result/N20230320ProteomEx_v2_demo_whole_library_Yueliang/Whole_lib_report.gg_matrix.tsv')
matgg <- dfgg %>% column_to_rownames('Genes') %>% log2() # log2 transform
colnames(matgg) %<>% str_extract('b\\d+_(\\d+|pool)')
dfgg <- matgg %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .) %>% 
  filter(Batch %in% dfpg$Batch)
matgg <- matgg[, dfgg$Batch[!is.na(dfgg$Region)]]
dfgg <- matgg %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .)

# Missing ratio cut-off: <20% in at least 3 Regions
dfgg %>% count(Region)

dfgg_miss <- dfgg %>% group_by(Region) %>% select(-(SampleName:BatchHead)) %>%
  summarise_all(function(y) sum(is.na(y)) / length(y)) %>% 
  mutate(Region = as.character(Region))
dfgg_miss %<>% rbind(c(Region = 'ALL', apply(matgg, 1, function(x) sum(is.na(x)) / length(x))))

tmp <- dfgg_miss %>%
  filter(Region != 'ALL') %>% 
  column_to_rownames('Region') %>%
  apply(2, function(y) sum(y < 0.2) >= 3 )
selected_gg <- names(tmp[tmp])
dfgg_miss %<>% select(Region, all_of(selected_gg))
matgg <- matgg[colnames(dfgg_miss)[-1], ]

### missing ratio control -----
mrgg <- apply(matgg, 1, function(x) sum(is.na(x)) / length(x))
plot(density(mrgg))

x <- seq(0.01, 0.99, 0.01)
y <- sapply(x, function(mrgg_cutoff){
  X <- matgg[mrgg < mrgg_cutoff, ]
  return(sum(is.na(X)) / nrow(X) / ncol(X))
})
data.frame(x, y)
plot(x * 100, y * 100, xlab = '% NA ratio cutoff on protein level', ylab = '% whole matggrix NA ratio')
abline(h = 6.5, v = 50)

# dfgg_miss %>% select(Region, all_of(colnames(dfgg_miss)[-1][dfgg_miss[7, -1] < 0.5])) %>% identical(dfgg_miss) # TRUE

matgg <- matgg[colnames(dfgg_miss)[-1], ]
dfgg %<>% select(SampleName:BatchHead, all_of(rownames(matgg)))
dim(matgg) # 4101 123



library(GSVA)
library(GSEABase)
library(enrichplot)
library(limma)
library(pheatmap)
# H: hallmark gene sets (50 gene sets)
# C2_CP: Canonical pathways (3795 gene sets)
# C5_GO_BP: subset of Gene Ontology gene sets (7647 gene sets)
# C6: oncogenic signature gene sets (189 gene sets)
# C7: immunologic signature gene sets (5219 gene sets)
gsH <- getGmt("data/gmt/h.all.v2023.2.Hs.symbols.gmt")
gsC2CP <- getGmt("data/gmt/c2.cp.v2023.2.Hs.symbols.gmt")
gsC5GOBP <- getGmt("data/gmt/c5.go.bp.v2023.2.Hs.symbols.gmt")
gsC6 <- getGmt("data/gmt/c6.all.v2023.2.Hs.symbols.gmt")
gsC7 <- getGmt("data/gmt/c7.all.v2023.2.Hs.symbols.gmt")

### 1.7.1 C5_GO_BP: subset of Gene Ontology gene sets -----
GSVA_c5gobp <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gsC5GOBP,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)

identical(colnames(GSVA_c5gobp), dfgg$Batch) # TRUE

Type <- as.factor(dfgg$Region)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(
  contrasts = c('TypeL-TypeN', 'TypeH-TypeL', 'TypePC-TypeH', 'TypeCC-TypeH'),
  levels = design
)
fit1 <- lmFit(GSVA_c5gobp, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
# res_c5gobp <- topTable(fit3, number = Inf)
# qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
# abline(0,1)

res_c5gobp <- list(
  L_N = topTable(fit3, coef = 'TypeL-TypeN', number = Inf) %>% rownames_to_column('ID'),
  H_L = topTable(fit3, coef = 'TypeH-TypeL', number = Inf) %>% rownames_to_column('ID'),
  PC_H = topTable(fit3, coef = 'TypePC-TypeH', number = Inf) %>% rownames_to_column('ID'),
  CC_H = topTable(fit3, coef = 'TypeCC-TypeH', number = Inf) %>% rownames_to_column('ID')
) %>%
  plyr::ldply(.id = 'Coef')

res_c5gobp <- topTable(fit3, number = Inf) %>%
  dplyr::select(AveExpr:adj.P.Val) %>% 
  set_colnames(str_c('Total_', colnames(.))) %>% 
  rownames_to_column('ID') %>% 
  inner_join(res_c5gobp, .)

depw_c5gobp <- res_c5gobp %>% 
  dplyr::filter(Total_adj.P.Val < 0.05) %>%
  plyr::daply('ID', function(dfsub){
  ret <- dfsub %>% filter(adj.P.Val < 0.05)
  ifelse((all(ret$logFC > 0) | all(ret$logFC < 0)) &
           any(ret$logFC < quantile(res_c5gobp$logFC, 0.01) |
                 ret$logFC > quantile(res_c5gobp$logFC, 0.99)) &
           nrow(ret) > 1,
         T, F)
}) %>% which() %>% names()
res_c5gobp %>% filter(ID %in%depw_c5gobp) %>% dim()
res_c5gobp %>% filter(ID %in%depw_c5gobp) %>% count(ID) %>% dim()
x <- res_c5gobp %>% filter(ID %in% depw_c5gobp) %>% pull(ID) %>% unique()
geneIds(gsC5GOBP)[x] %>% sapply(length) %>% quantile()
str_replace_all(x, '_', ' ') %>% writeClipboard()
# str_subset(unique(res_c5gobp$ID), 'EXTRACELLULAR') %>% .[. %in% x]
# res_c5gobp[res_c5gobp$ID == 'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]
# GSVA_c5gobp['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]


X1 <- GSVA_c5gobp[x, dfheat$Batch]
rownames(X1) %<>% str_replace_all('_', ' ')
clust_col <- hclust(dist(t(X1), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
heatmap_GSVA_C5GOBP <- pheatmap(X1, scale = 'row',
         annotation_col = dfheat[, 'Region', drop = F],
         annotation_colors = anno_colors,
         cluster_rows = T, cluster_cols = clust_col,
         show_rownames = T, show_colnames = F,
         color = heat_colors, breaks = bk,
         fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
         border_color = F, na_col = '#AAAAAA',
         # gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
         filename = 'F5A_GSVA_C5GOBP_heatmap.pdf',
         width = 9, height = 4
)

### 1.7.2 H: hallmark gene sets -----
GSVA_h <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gsH,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)
identical(colnames(GSVA_h), dfgg$Batch) # TRUE

Type <- as.factor(dfgg$Region)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(
  contrasts = c('TypeL-TypeN', 'TypeH-TypeL', 'TypePC-TypeH', 'TypeCC-TypeH'),
  levels = design
)
fit1 <- lmFit(GSVA_h, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
# res_h <- topTable(fit3, number = Inf)
# qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
# abline(0,1)

res_h <- list(
  L_N = topTable(fit3, coef = 'TypeL-TypeN', number = Inf) %>% rownames_to_column('ID'),
  H_L = topTable(fit3, coef = 'TypeH-TypeL', number = Inf) %>% rownames_to_column('ID'),
  PC_H = topTable(fit3, coef = 'TypePC-TypeH', number = Inf) %>% rownames_to_column('ID'),
  CC_H = topTable(fit3, coef = 'TypeCC-TypeH', number = Inf) %>% rownames_to_column('ID')
) %>%
  plyr::ldply(.id = 'Coef')

res_h <- topTable(fit3, number = Inf) %>%
  dplyr::select(AveExpr:adj.P.Val) %>% 
  set_colnames(str_c('Total_', colnames(.))) %>% 
  rownames_to_column('ID') %>% 
  inner_join(res_h, .)

depw_h <- res_h %>% 
  dplyr::filter(Total_adj.P.Val < 0.05) %>%
  plyr::daply('ID', function(dfsub){
    ret <- dfsub %>% filter(adj.P.Val < 0.05)
    ifelse((all(ret$logFC > 0) | all(ret$logFC < 0)) &
             # any(ret$logFC < quantile(res_h$logFC, 0.01) |
             #       ret$logFC > quantile(res_h$logFC, 0.99)) &
             nrow(ret) > 1,
           T, F)
  }) %>% which() %>% names()
res_h %>% filter(ID %in%depw_h) %>% dim()
res_h %>% filter(ID %in%depw_h) %>% count(ID) %>% dim()
x <- res_h %>% filter(ID %in% depw_h) %>% pull(ID) %>% unique()
geneIds(gsH)[x] %>% sapply(length) %>% quantile()
str_replace_all(x, '_', ' ') %>% writeClipboard()
# str_subset(unique(res_h$ID), 'EXTRACELLULAR') %>% .[. %in% x]
# res_h[res_h$ID == 'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]
# GSVA_h['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]

res_h %>% filter(ID %in% x) %>% arrange(ID)
X2 <- GSVA_h[x, dfheat$Batch]
rownames(X2) %<>% str_replace_all('_', ' ')
heatmap_GSVA_H <- pheatmap(X2, scale = 'row',
                         annotation_col = dfheat[, 'Region', drop = F],
                         annotation_colors = anno_colors,
                         cluster_rows = T, cluster_cols = F,
                         show_rownames = T, show_colnames = F,
                         color = heat_colors, breaks = bk,
                         fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                         border_color = F, na_col = '#AAAAAA',
                         gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                         filename = 'F5A_GSVA_Hallmark_heatmap.pdf',
                         width = 9, height = 0.8
)

### 1.7.3 C2_CP: Canonical pathways -----
GSVA_c2cp <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gsC2CP,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)

identical(colnames(GSVA_c2cp), dfgg$Batch) # TRUE

Type <- as.factor(dfgg$Region)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(
  contrasts = c('TypeL-TypeN', 'TypeH-TypeL', 'TypePC-TypeH', 'TypeCC-TypeH'),
  levels = design
)
fit1 <- lmFit(GSVA_c2cp, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
# res_c2cp <- topTable(fit3, number = Inf)
# qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
# abline(0,1)

res_c2cp <- list(
  L_N = topTable(fit3, coef = 'TypeL-TypeN', number = Inf) %>% rownames_to_column('ID'),
  H_L = topTable(fit3, coef = 'TypeH-TypeL', number = Inf) %>% rownames_to_column('ID'),
  PC_H = topTable(fit3, coef = 'TypePC-TypeH', number = Inf) %>% rownames_to_column('ID'),
  CC_H = topTable(fit3, coef = 'TypeCC-TypeH', number = Inf) %>% rownames_to_column('ID')
) %>%
  plyr::ldply(.id = 'Coef')

res_c2cp <- topTable(fit3, number = Inf) %>%
  dplyr::select(AveExpr:adj.P.Val) %>% 
  set_colnames(str_c('Total_', colnames(.))) %>% 
  rownames_to_column('ID') %>% 
  inner_join(res_c2cp, .)

depw_c2cp <- res_c2cp %>% 
  dplyr::filter(Total_adj.P.Val < 0.01) %>%
  plyr::daply('ID', function(dfsub){
    ret <- dfsub %>% filter(adj.P.Val < 0.01)
    ifelse((all(ret$logFC > 0) | all(ret$logFC < 0)) &
             any(ret$logFC < quantile(res_c2cp$logFC, 0.05) |
                   ret$logFC > quantile(res_c2cp$logFC, 0.95)) &
             nrow(ret) > 1,
           T, F)
  }) %>% which() %>% names()
res_c2cp %>% filter(ID %in%depw_c2cp) %>% dim()
res_c2cp %>% filter(ID %in%depw_c2cp) %>% count(ID) %>% dim()
x <- res_c2cp %>% filter(ID %in% depw_c2cp) %>% pull(ID) %>% unique()
geneIds(gsC2CP)[x] %>% sapply(length) %>% quantile()
str_replace_all(x, '_', ' ') %>% writeClipboard()
# str_subset(unique(res_c2cp$ID), 'EXTRACELLULAR') %>% .[. %in% x]
# res_c2cp[res_c2cp$ID == 'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]
# GSVA_c2cp['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]


X3 <- GSVA_c2cp[x, dfheat$Batch, drop = F]
rownames(X3) %<>% str_replace_all('_', ' ')
# clust_col <- hclust(dist(t(X3), method = 'euclidean'), method = 'ward.D2') %>% 
#   as.dendrogram() %>% 
#   reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
#   as.hclust()
heatmap_GSVA_C2CP <- pheatmap(X3, scale = 'row',
                              annotation_col = dfheat[, 'Region', drop = F],
                              annotation_colors = anno_colors,
                              cluster_rows = F, cluster_cols = F,
                              show_rownames = T, show_colnames = F,
                              color = heat_colors, breaks = bk,
                              fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                              border_color = F, na_col = '#AAAAAA',
                              gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                              filename = 'F5A_GSVA_C2CP_heatmap.pdf',
                              width = 20, height = 1.2
)

### 1.7.4 C6: oncogenic signature gene sets -----
GSVA_c6 <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gsC6,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)

identical(colnames(GSVA_c6), dfgg$Batch) # TRUE

Type <- as.factor(dfgg$Region)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(
  contrasts = c('TypeL-TypeN', 'TypeH-TypeL', 'TypePC-TypeH', 'TypeCC-TypeH'),
  levels = design
)
fit1 <- lmFit(GSVA_c6, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
# res_c6 <- topTable(fit3, number = Inf)
# qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
# abline(0,1)

res_c6 <- list(
  L_N = topTable(fit3, coef = 'TypeL-TypeN', number = Inf) %>% rownames_to_column('ID'),
  H_L = topTable(fit3, coef = 'TypeH-TypeL', number = Inf) %>% rownames_to_column('ID'),
  PC_H = topTable(fit3, coef = 'TypePC-TypeH', number = Inf) %>% rownames_to_column('ID'),
  CC_H = topTable(fit3, coef = 'TypeCC-TypeH', number = Inf) %>% rownames_to_column('ID')
) %>%
  plyr::ldply(.id = 'Coef')

res_c6 <- topTable(fit3, number = Inf) %>%
  dplyr::select(AveExpr:adj.P.Val) %>% 
  set_colnames(str_c('Total_', colnames(.))) %>% 
  rownames_to_column('ID') %>% 
  inner_join(res_c6, .)

depw_c6 <- res_c6 %>% 
  dplyr::filter(Total_adj.P.Val < 0.05) %>%
  plyr::daply('ID', function(dfsub){
    ret <- dfsub# %>% filter(adj.P.Val < 0.1)
    ifelse((all(ret$logFC > 0) | all(ret$logFC < 0)) &
             any(abs(ret$logFC) > 0.05) &
             # any(ret$logFC < quantile(res_c6$logFC, 0.01) |
             #       ret$logFC > quantile(res_c6$logFC, 0.99)) &
             nrow(ret) > 1,
           T, F)
  }) %>% which() %>% names()
res_c6 %>% filter(ID %in%depw_c6) %>% dim()
res_c6 %>% filter(ID %in%depw_c6) %>% count(ID) %>% dim()
x <- res_c6 %>% filter(ID %in% depw_c6) %>% pull(ID) %>% unique()
geneIds(gsC6)[x] %>% sapply(length) %>% quantile()
str_replace_all(x, '_', ' ') %>% writeClipboard()
# str_subset(unique(res_c6$ID), 'EXTRACELLULAR') %>% .[. %in% x]
# res_c6[res_c6$ID == 'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]
# GSVA_c6['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]


bk <- unique(c(seq(-2, 2, length = 50)))
X4 <- GSVA_c6[x, dfheat$Batch, drop = F]
rownames(X4) %<>% str_replace_all('_', ' ')
heatmap_GSVA_C6 <- pheatmap(X4, scale = 'row',
                            annotation_col = dfheat[, 'Region', drop = F],
                            annotation_colors = anno_colors,
                            cluster_rows = F, cluster_cols = F,
                            show_rownames = T, show_colnames = F,
                            color = heat_colors, breaks = bk,
                            fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                            border_color = F, na_col = '#AAAAAA',
                            # gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                            filename = 'F5A_GSVA_C6_heatmap.pdf',
                            width = 9, height = 0.6
)

### 1.7.5 C7: immunologic signature gene sets -----
GSVA_c7 <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gsC7,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)

identical(colnames(GSVA_c7), dfgg$Batch) # TRUE

Type <- as.factor(dfgg$Region)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(
  contrasts = c('TypeL-TypeN', 'TypeH-TypeL', 'TypePC-TypeH', 'TypeCC-TypeH'),
  levels = design
)
fit1 <- lmFit(GSVA_c7, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
# res_c7 <- topTable(fit3, number = Inf)
# qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
# abline(0,1)

res_c7 <- list(
  L_N = topTable(fit3, coef = 'TypeL-TypeN', number = Inf) %>% rownames_to_column('ID'),
  H_L = topTable(fit3, coef = 'TypeH-TypeL', number = Inf) %>% rownames_to_column('ID'),
  PC_H = topTable(fit3, coef = 'TypePC-TypeH', number = Inf) %>% rownames_to_column('ID'),
  CC_H = topTable(fit3, coef = 'TypeCC-TypeH', number = Inf) %>% rownames_to_column('ID')
) %>%
  plyr::ldply(.id = 'Coef')

res_c7 <- topTable(fit3, number = Inf) %>%
  dplyr::select(AveExpr:adj.P.Val) %>% 
  set_colnames(str_c('Total_', colnames(.))) %>% 
  rownames_to_column('ID') %>% 
  inner_join(res_c7, .)

depw_c7 <- res_c7 %>% 
  dplyr::filter(Total_adj.P.Val < 0.05) %>%
  plyr::daply('ID', function(dfsub){
    ret <- dfsub %>% filter(adj.P.Val < 0.05)
    ifelse((all(ret$logFC > 0) | all(ret$logFC < 0)) &
             any(ret$logFC < quantile(res_c7$logFC, 0.01) |
                   ret$logFC > quantile(res_c7$logFC, 0.99)) &
             nrow(ret) > 1,
           T, F)
  }) %>% which() %>% names()
res_c7 %>% filter(ID %in%depw_c7) %>% dim()
res_c7 %>% filter(ID %in%depw_c7) %>% count(ID) %>% dim()
x <- res_c7 %>% filter(ID %in% depw_c7) %>% pull(ID) %>% unique()
geneIds(gsC7)[x] %>% sapply(length) %>% quantile()
str_replace_all(x, '_', ' ') %>% writeClipboard()
# str_subset(unique(res_c7$ID), 'EXTRACELLULAR') %>% .[. %in% x]
# res_c7[res_c7$ID == 'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]
# GSVA_c7['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]


X5 <- GSVA_c7[x, dfheat$Batch]
rownames(X5) %<>% str_replace_all('_', ' ')

clust_col <- hclust(dist(t(X5), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
heatmap_GSVA_C7 <- pheatmap(X5, scale = 'row',
                            annotation_col = dfheat[, 'Region', drop = F],
                            annotation_colors = anno_colors,
                            cluster_rows = T, cluster_cols = F,#cluster_cols = clust_col,
                            show_rownames = T, show_colnames = F,
                            color = heat_colors, breaks = bk,
                            fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                            border_color = F, na_col = '#AAAAAA',
                            gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                            filename = 'F5A_GSVA_C7_heatmap.pdf',
                            width = 12, height = 6
)


### 1.7.6 LM22 -----
load('data/gmt/gene_set_LM22.Rdata')

GSVA_lm22 <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gs_LM22,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)
unlist(gs_LM22) %in% rownames(matgg) %>% sum() # 42
unlist(gs_LM22) %in% rownames(matgg) %>% length() # 1140



## 1.8 GSEA --------------------
# prepare gene list
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

# use all detected genes
DEA_gsea <- DEA %>% 
  filter(Comparison %in% my_comparisons) %>%
  drop_na(P, Log2FC) %>%
  arrange(P) %>%
  plyr::ddply('Protein', function(dfsub){
    ret <- dfsub %>% 
      slice(1) %>% select(Protein, Genes, Log2FC, P)
    ret$Consistent <- all(dfsub$Log2FC > 0) | all(dfsub$Log2FC < 0)
    return(ret)
  })
DEA_gsea %<>% arrange(desc(Log2FC))
# DEA_gsea_sig <- DEA_gsea %>% filter(P < 0.05, Consistent)

tmp1 <- clusterProfiler::bitr(DEA_gsea$Genes, OrgDb = 'org.Hs.eg.db', fromType = 'SYMBOL', toType = 'ENTREZID')
tmp2 <- clusterProfiler::bitr(DEA_gsea$Protein, OrgDb = 'org.Hs.eg.db', fromType = 'UNIPROT', toType = 'ENTREZID')
tmp <- tmp1 %>% full_join(tmp2) %>% setNames(c('Genes', 'ENTREZID', 'Protein'))
DEA_gsea %<>% left_join(tmp)

gene_list <- drop_na(DEA_gsea, Genes)$Log2FC %>% setNames(drop_na(DEA_gsea, Genes)$Genes)
gid_list <- drop_na(DEA_gsea, ENTREZID)$Log2FC %>% setNames(drop_na(DEA_gsea, ENTREZID)$ENTREZID)


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

### GOBP -----
gseGOBP <- gseGO(gene_list, OrgDb = org.Hs.eg.db,
                 ont = 'BP', keyType = 'SYMBOL',
                 minGSSize = 10, maxGSSize = 3000,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
df_gsegobp <- gseGOBP[]

# goplot(gseGOBP, showCategory = 10)
# gseaplot2(gseGOBP, 'GO:0090304', pvalue_table = T, ES_geom = 'line')
plot_gseGOBP <- gseGOBP %>% filter(p.adjust < 1e-3, abs(NES) > 1.8, setSize < 300) %>% ridgeplot(showCategory = 10, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5A_DEPs_progress_109_GSEA_GOBP.pdf', plot_gseGOBP, width = 8, height = 6)

### GOCC -----
gseGOCC <- gseGO(gene_list, OrgDb = org.Hs.eg.db,
                 ont = 'CC', keyType = 'SYMBOL',
                 minGSSize = 30, maxGSSize = 5000,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
df_gsegocc <- gseGOCC[]

# goplot(gseGOCC, showCategory = 10)
# gseaplot2(gseGOCC, 'GO:0090304', pvalue_table = T, ES_geom = 'line')
plot_gseGOCC <- gseGOCC %>% filter(p.adjust < 1e-3, abs(NES) > 1.8) %>% ridgeplot(showCategory = 10, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5A_DEPs_progress_109_GSEA_GOCC.pdf', plot_gseGOCC, width = 8, height = 6)

### GOMF -----
gseGOMF <- gseGO(gene_list, OrgDb = org.Hs.eg.db,
                 ont = 'MF', keyType = 'SYMBOL',
                 minGSSize = 30, maxGSSize = 5000,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
df_gsegomf <- gseGOMF[]

# goplot(gseGOMF, showCategory = 10)
# gseaplot2(gseGOMF, 'GO:0090304', pvalue_table = T, ES_geom = 'line')
plot_gseGOMF <- gseGOMF %>%
  # filter(p.adjust < 1e-3, abs(NES) > 1.8) %>%
  ridgeplot(showCategory = 10, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5A_DEPs_progress_109_GSEA_GOMF.pdf', plot_gseGOMF, width = 8, height = 6)

### KEGG -----
library(KEGG.db)

gseKEGG <- gseKEGG(gid_list, organism = 'hsa', use_internal_data = T,
                   keyType = 'ncbi-geneid',
                   minGSSize = 5, maxGSSize = 2500,
                   exponent = 1, eps = 0,
                   pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                   verbose = T, seed = F, by = 'fgsea')
df_gsekegg <- gseKEGG[]

library(KEGGREST)
Classes <- lapply(df_gsekegg$ID, function(KEGGID){
  keggGet(KEGGID)[[1]]$CLASS
})
KEGGID_filtered <- lapply(Classes, function(x) str_detect(x, 'Human Diseases') & !str_detect(x, '[Cc]ancer')) %>% unlist() %>% .[!.] %>% names()
df_gsekegg %>% filter(ID %in% KEGGID_filtered)

# gseaplot2(gseKEGG, 'hsa04151', pvalue_table = T, ES_geom = 'line')
# gseaplot2(gseKEGG, 'hsa00190', pvalue_table = T, ES_geom = 'line')
# gseaplot2(gseKEGG, 'hsa04512', pvalue_table = T, ES_geom = 'line')
plot_gseKEGG <- gseKEGG %>%
  filter(ID %in% KEGGID_filtered) %>% 
  # filter(p.adjust < 1e-3, abs(NES) > 1.8) %>%
  ridgeplot(#showCategory = 10,
            orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5A_GSEA_progress_KEGG.pdf', plot_gseKEGG, width = 8, height = 8)

pdf('F5A_GSEA_progress_KEGG_PI3K-Akt.pdf', width = 7, height = 7)
gseaplot2(gseKEGG, 'hsa04151', title = 'KEGG: PI3K-Akt signaling pathway', pvalue_table = T, ES_geom = 'line')
graphics.off()


## 1.9 Protein expression -----
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

### 1.9.1 KEGG; semantic similarity analysis ------
#### hsa04151 PI3K-Akt signaling pathway -----
##### heatmap ---------
gid_hsa04151 <- dfkegg %>% filter(ID == 'hsa04151') %>% pull(geneID) %>%
  str_split('/') %>% .[[1]]

df_hsa04151 <- clusterProfiler::bitr(gid_hsa04151, OrgDb = 'org.Hs.eg.db', fromType = 'ENTREZID', toType = 'SYMBOL') %>% 
  rename(GeneID = ENTREZID, Genes = SYMBOL) %>% 
  left_join(DEP_prog, )
prot_hsa04151 <- unique(df_hsa04151$Protein)

library(pheatmap)
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

clust_col <- hclust(dist(t(mat_heat[prot_hsa04151, ]), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
heat_hsa04151 <- pheatmap(mat_heat[prot_hsa04151, ],
                          annotation_col = dfheat %>% select(Region),
                          annotation_colors = anno_colors,
                          cluster_rows = T, cluster_cols = clust_col,
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
X <- mat_heat[heat_hsa04151$tree_row$labels[heat_hsa04151$tree_row$order], ]
rownames(X) <- DEP_prog %>% distinct(Protein, Label) %>% set_rownames(.$Protein) %>% .[rownames(X), ] %>% pull(Label)
pF_heat3 <- pheatmap(X,
                     annotation_col = dfheat %>% select(Region),
                     annotation_colors = anno_colors,
                     cluster_rows = F, cluster_cols = F,
                     show_rownames = T, show_colnames = F,
                     color = heat_colors, breaks = bk,
                     fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                     border_color = F, na_col = '#AAAAAA',
                     gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                     filename = 'F5A_DEPs_progress_KEGGhsa04151_PI3K_Akt_heatmap.pdf',
                     width = 4, height = 1.5
)

##### boxplot ---------
dfbox <- dfpg %>%
  filter(!is.na(Region)) %>%
  select(SampleName:BatchHead, all_of(prot_hsa04151)) %>%
  pivot_longer(cols = -c(SampleName:BatchHead), names_to = 'Protein', values_to = 'Log2Intensity') %>%
  left_join(DEA %>% distinct(Protein, Genes, Label), relationship = 'many-to-many') %>%
  mutate(Protein = factor(Protein, levels = prot_hsa04151))

boxplot_ls <- plyr::dlply(dfbox, 'Protein', function(dfsub){
  label <- dfsub$Label[1]
  set.seed(10)
  ggplot(dfsub, aes(x = Region, y = `Log2Intensity`, color = Region)) +
    stat_boxplot(geom = 'errorbar', width = 0.4)+
    geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
    geom_jitter(alpha = 0.4, size = 1, width = 0.15)+
    labs(x = '', y = '', subtitle = label) +
    scale_color_manual(values = mycolors) +
    theme_classic() +
    theme(text = element_text(size = 15), legend.position = 'none') +
    ggpubr::stat_compare_means(
      method = 't.test',
      map_signif_level = F,# hide.ns = F, label = 'p.signif',
      comparisons = str_split(my_comparisons, ' vs. '),
      size = 5, hjust = 0.5, vjust = 0)
})
boxplots <- ggpubr::ggarrange(plotlist = boxplot_ls, nrow = 3, ncol = 3)
ggsave('F5A_DEPs_progress_KEGGhsa04151_PI3K_Akt_box.pdf', boxplots, width = 3 * 4, height = 3 * 4, limitsize = F)

list(heatmap = X %>% as.data.frame() %>% rownames_to_column('Label'),
     annotation_col = dfheat %>% select(Region),
     box = dfbox) %>% 
rio::export('F5A_DEPs_progress_KEGGhsa04151_PI3K_Akt_box.xlsx')



## 1.x output data --------------------
list(#DEP_prog109 = DEP_prog,
     DEP_prog50 = DEP_prog2,
     DEP_prog50_expr = dfbox,
     DEP_prog50_heatmap = X %>% as.data.frame() %>% rownames_to_column(),
     DEP_prog50_UMAP = df_umap) %>%
  rio::export('F5A_DEP_disease_progress.xlsx')

pkgs <- data.frame(Package = c('clusterProfiler', 'ReactomePA', 'org.Hs.eg.db'))
pkgs$Version <- sapply(pkgs$Package, function(pkg) as.character(packageVersion(pkg)))
list(#GOBP_sim_sim = GOBP_simplify_sim %>% arrange(cluster) %>% inner_join(dfgobp_simplify, ., by = c(ID = 'id')),
     # GOBP_similar_mat = similar_matrix %>% as.data.frame() %>% rownames_to_column(),
     # KEGG = dfkegg,
     Reactome = dfreactome,
     PackageVersion = pkgs) %>%
  rio::export('F5A_DEPs_progress_50_pathway_enrichment.xlsx')

# pkgs <- data.frame(Package = c('GSVA', 'GSEABase', 'limma'))
# pkgs$Version <- sapply(pkgs$Package, function(pkg) as.character(packageVersion(pkg)))
# list(GSVA_c5gobp = GSVA_c5gobp %>% as.data.frame() %>% rownames_to_column('ID'),
#      GSVA_h = GSVA_h %>% as.data.frame() %>% rownames_to_column('ID'),
#      GSVA_c2cp = GSVA_c2cp %>% as.data.frame() %>% rownames_to_column('ID'),
#      GSVA_c6 = GSVA_c6 %>% as.data.frame() %>% rownames_to_column('ID'),
#      GSVA_c7 = GSVA_c7 %>% as.data.frame() %>% rownames_to_column('ID'),
#      res_c5gobp = res_c5gobp,
#      res_h = res_h,
#      res_c2cp = res_c2cp,
#      res_c6 = res_c6,
#      res_c7 = res_c7,
#      pkgs = pkgs) %>% 
#   rio::export('F5A_progress_GSVA.xlsx')
# 
# 
# list(GSEA_GOBP = df_gsegobp,
#      GSEA_GOCC =  df_gsegocc,
#      GSEA_GOMF = df_gsegomf,
#      GSEA_KEGG = df_gsekegg) %>% 
#   rio::export('F5A_progress_GSEA.xlsx')

save.image('Figure5A_20240107.RData')



# 2.Others vs. N ----------------------------------------------------------
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
rm(list = ls()); gc()
load('../Figure4/Figure4.RData') # based on Figure4 data

## 2.1 statistics ------------------
### figure4 results ----
dim(DEA) # 67510 16
DEP %>% pull(Protein) %>% unique() %>% length() # 555 DEPs in total
DEP_solid %>% pull(Protein) %>% unique() %>% length() # 516


### cancer associated proteins -----
my_comparisons <- c('L vs. N', 'H vs. N', 'C vs. N', 'PC vs. N', 'CC vs. N')
DEA_cancer <- DEA %>% filter(Comparison %in% my_comparisons,
                             adj.P < 0.05)
prot_cancer <- DEA_cancer %>%
  pivot_wider(id_cols = 'Protein', names_from = 'Comparison', values_from = 'Log2FC') %>% 
  column_to_rownames('Protein') %>% 
  mutate_all(function(x) ifelse(is.na(x), 0, x)) %>% # do not consider NA fold-change
  apply(1, function(x) (all(x >= 0) | all(x <= 0)) & any(abs(x) > log2(2))) %>% # fold-change values should be in the same trend, and at least one of them >2
  which() %>% names()
DEP_cancer <- DEA_cancer %>% filter(Protein %in% prot_cancer)


### boxplot ------
prot_cancer1 <- DEP_cancer %>% count(Protein) %>% filter(n == length(my_comparisons)) %>% pull(Protein) # all five significant comparison
DEP_cancer1 <- DEP_cancer %>% filter(Protein %in% prot_cancer1)

prot_cancer2 <- DEP_cancer1 %>% select(matches('^Missing.+_')) %>% apply(1, function(x) all(x < 0.5)) %>% filter(DEP_cancer1, .) %>% pull(Protein) %>% unique() # Missing <50% in all Regions
DEP_cancer2 <- DEP_cancer1 %>% filter(Protein %in% prot_cancer2)

dfbox <- dfpg %>%
  filter(!is.na(Region)#, Region != 'C'
         ) %>% 
  select(SampleName:BatchHead, all_of(prot_cancer2)) %>% 
  pivot_longer(cols = -c(SampleName:BatchHead), names_to = 'Protein', values_to = 'Log2Intensity') %>% 
  left_join(DEA %>% distinct(Protein, Genes, Label), relationship = 'many-to-many') %>%
  mutate(Protein = factor(Protein, levels = prot_cancer2))

nafill <- dfpg %>% select(-(SampleName:BatchHead)) %>%
  min(na.rm = T) %>% `+`(log2(0.8))
boxplot_ls <- plyr::dlply(dfbox, 'Protein', function(dfsub){
  dfna <- dfsub %>% filter(is.na(Log2Intensity)) %>% 
    mutate(Log2Intensity = nafill)
  ymax <- max(dfsub$Log2Intensity, na.rm = T)
  DEPsub <- DEA %>%
    filter(
      Protein == dfsub$Protein[1], adj.P < 0.05,
      Comparison %in% my_comparisons
    ) %>%
    mutate(y = ymax * seq(1.05, by = 0.08, length.out = nrow(.)))
  label <- dfsub$Label[1]
  set.seed(10)
  ggplot(dfsub, aes(x = Region, y = `Log2Intensity`, color = Region)) +
    stat_boxplot(geom = 'errorbar', width = 0.4)+
    geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
    geom_jitter(alpha = 0.4, size = 1, width = 0.15)+
    geom_jitter(data = dfna, alpha = 0.4, size = 1, width = 0.15) +
    labs(x = '', y = '', subtitle = label) +
    scale_y_continuous(limits = c(0, ymax * 1.3)) +
    scale_color_manual(values = mycolors) +
    theme_classic() +
    theme(text = element_text(size = 15), legend.position = 'none') +
    # ggpubr::stat_compare_means(
    #   method = 't.test',
    #   map_signif_level = F, hide.ns = F, label = 'p.signif',
    #   comparisons = str_split(my_comparisons, ' vs. '),
    #   size = 5, hjust = 0.5, vjust = 0) +
    ggsignif::geom_signif(
      data = DEPsub,
      manual = T, inherit.aes = F,
      aes(xmax = str_remove(Comparison, '^.+vs\\. '),
          xmin = str_remove(Comparison, ' vs\\. .+$'),
          y_position = y,
          annotations = sprintf('%.2e', adj.P))
    )
})
boxplots <- ggpubr::ggarrange(plotlist = boxplot_ls, nrow = 41, ncol = 8)
ggsave('F5B_DEPs_cancerous_box.pdf', boxplots, width = 8 * 4, height = 41 * 4, limitsize = F)

## 2.2 volcano ------------
set.seed(10)
p5B_volcano <- DEA %>% filter(Comparison %in% my_comparisons) %>% 
  ggplot(aes(x = Comparison, y = Log2FC, color = interaction(Type, is.significant)))+
  geom_jitter(#data = DEA %>% filter(abs(Log2FC) > 2),
    alpha = 0.5, size = 0.5, width = 0.4) +
  labs(x = '', y = 'Log2 Fold-Change', subtitle = '# differentially expressed proteins')+
  scale_color_manual(name = '',
                     values = c(Up.FALSE = '#DDDDDD', Down.FALSE = '#DDDDDD',
                                Up.TRUE = 'red3', Down.TRUE = 'blue3'),
                     labels = c('', '|FC|<2 or adj.P>0.05', 'Down', 'Up')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'top') +
  # geom_jitter(data = DEP_cancer2 %>% filter(is.significant),
  #             size = 3, width = 0.4) +
  ggrepel::geom_text_repel(
    data = DEP_cancer2 %>% filter(is.significant),
    aes(label = Genes, size = abs(Log2FC)),
    direction = 'x', show.legend = F,
    force = 1.2, seed = 10, max.overlaps = 100, segment.color = NA
  ) +
  scale_size_continuous(range = c(0.8, 4), breaks = c(0.5, 3))
ggsave('F5B_DEPs_cancerous_volcano.pdf', p5B_volcano, width = 14, height = 8)


## 2.3 heatmap ------
library(pheatmap)
mat_heat <- dfheat %>%
  select(-(Batch:BatchHead)) %>% 
  scale() %>% t() # Z-score
mat_heat[is.na(mat_heat)] <- min(mat_heat, na.rm = T)

clust_col <- hclust(dist(t(mat_heat[unique(DEP_cancer2$Protein), ]), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
heatTN <- pheatmap(mat_heat[unique(DEP_cancer2$Protein), ],
                   annotation_col = dfheat %>% select(Region),
                   annotation_colors = anno_colors,
                   cluster_rows = T, cluster_cols = clust_col,
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
X <- mat_heat[heatTN$tree_row$labels[heatTN$tree_row$order], ]
rownames(X) <- DEP_cancer2 %>% distinct(Protein, Label) %>% set_rownames(.$Protein) %>% .[rownames(X), ] %>% pull(Label)
p5B_heat <- pheatmap(X,
                     annotation_col = dfheat %>% select(Region),
                     annotation_colors = anno_colors,
                     cluster_rows = F, cluster_cols = clust_col,
                     show_rownames = T, show_colnames = F,
                     color = heat_colors, breaks = bk,
                     fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                     border_color = F, na_col = '#AAAAAA',
                     # gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                     filename = 'F5B_DEPs_cancerous_heatmap.pdf',
                     width = 10, height = 8
)

## 2.4 UMAP ---------
library(umap)
X_imp <- X
X_imp[is.na(X_imp)] <- 0

umap <- umap(t(X_imp), n_neighbors = 60, min_dist = 0.8) # default n_neighbors = 10, min_dist = 0.1
df_umap <- info %>% inner_join(umap$layout %>% as.data.frame() %>% rownames_to_column('Batch'))
df_umap_center <- df_umap %>% group_by(Region) %>% summarise_at(vars(center_x = V1, center_y = V2), mean) # calculate center of clusters
df_umap %<>% left_join(df_umap_center)

p_umap <- ggplot(df_umap, aes(x = V1, y = V2, color = Region)) +
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
ggsave('F5B_DEPs_cancerous_UMAP.pdf', p_umap, width = 4, height = 4)

p_umap <- ggplot(df_umap, aes(x = V1, y = V2, color = Patient)) +
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
  theme(text = element_text(size = 10))
ggsave('F5B_DEPs_cancerous_UMAP_patient.pdf', p_umap, width = 4.2, height = 4)


## 2.5 Pathway enrichment -----------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

### Metascape -------
writeClipboard(unique(DEP_cancer2$Genes))
up_gene <- DEP_cancer2$Genes[DEP_cancer2$Type == 'Up'] %>% unique() # 132
down_gene <- DEP_cancer2$Genes[DEP_cancer2$Type == 'Down'] %>% unique() # 195


### GO ------------
#### Gene ontology profile ----
GOup <- enrichGO(unique(up_gene), OrgDb = org.Hs.eg.db,
                 ont = 'ALL', keyType = 'SYMBOL', readable = T,
                 minGSSize = 5, maxGSSize = 5000, pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOup_simplify <- clusterProfiler::simplify(GOup, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgoup <- GOup[]
dfgoup_simplify <- GOup_simplify[]

GOdown <- enrichGO(unique(down_gene), OrgDb = org.Hs.eg.db,
                   ont = 'ALL', keyType = 'SYMBOL', readable = T,
                   minGSSize = 5, maxGSSize = 5000, pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOdown_simplify <- clusterProfiler::simplify(GOdown, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgodown <- GOdown[]
dfgodown_simplify <- GOdown_simplify[]

enrichplot::dotplot(GOup, showCategory = 10, split = 'ONTOLOGY', title = 'Gene Ontology') + facet_grid(ONTOLOGY~., scale='free')
enrichplot::dotplot(GOdown, showCategory = 10, split = 'ONTOLOGY', title = 'Gene Ontology') + facet_grid(ONTOLOGY~., scale='free')

#### focus on GOBP -----
GOBPup <- enrichGO(up_gene, OrgDb = org.Hs.eg.db,
                   ont = 'BP', keyType = 'SYMBOL', readable = T,
                   minGSSize = 50, maxGSSize = 2000, pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOBPup_simplify <- clusterProfiler::simplify(GOBPup, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgobpup <- GOBPup[]
dfgobpup_simplify <- GOBPup_simplify[]

GOBPdown <- enrichGO(down_gene, OrgDb = org.Hs.eg.db,
                     ont = 'BP', keyType = 'SYMBOL', readable = T,
                     minGSSize = 50, maxGSSize = 2000, pAdjustMethod = 'BH',
                     pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOBPdown_simplify <- clusterProfiler::simplify(GOBPdown, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgobpdown <- GOBPdown[]
dfgobpdown_simplify <- GOBPdown_simplify[]

plot_gobpup <- enrichplot::dotplot(GOBPup_simplify, showCategory = 8, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'GOBP (simplified) up') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_gobpdown <- enrichplot::dotplot(GOBPdown_simplify, showCategory = 8, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'GOBP (simplified) down') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_gobp <- ggpubr::ggarrange(plot_gobpup, plot_gobpdown, nrow = 1, ncol = 2)
ggsave('F5B_DEPs_cancerous_328_GOBP.pdf', plot_gobp, width = 8.5, height = 3.5)

pdf('F5B_DEPs_cancer_up_132_GOBP_similar.pdf', width = 10, height = 8)
similar_matrix_up <- simplifyEnrichment::GO_similarity(dfgobpup_simplify$ID, ont = 'BP', db = 'org.Hs.eg.db')
GOBPup_simplify_sim <- simplifyEnrichment::simplifyGO(similar_matrix_up)
graphics.off()

pdf('F5B_DEPs_cancer_down_195_GOBP_similar.pdf', width = 10, height = 8)
similar_matrix_down <- simplifyEnrichment::GO_similarity(dfgobpdown_simplify$ID, ont = 'BP', db = 'org.Hs.eg.db')
GOBPdown_simplify_sim <- simplifyEnrichment::simplifyGO(similar_matrix_down)
graphics.off()

#### focus on GOCC -----
GOCCup <- enrichGO(up_gene, OrgDb = org.Hs.eg.db,
                   ont = 'CC', keyType = 'SYMBOL', readable = T,
                   minGSSize = 50, maxGSSize = 5000, pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOCCup_simplify <- clusterProfiler::simplify(GOCCup, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgoccup <- GOCCup[]
dfgoccup_simplify <- GOCCup_simplify[]

GOCCdown <- enrichGO(down_gene, OrgDb = org.Hs.eg.db,
                     ont = 'CC', keyType = 'SYMBOL', readable = T,
                     minGSSize = 50, maxGSSize = 5000, pAdjustMethod = 'BH',
                     pvalueCutoff = 0.05, qvalueCutoff = 0.05)
GOCCdown_simplify <- clusterProfiler::simplify(GOCCdown, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgoccdown <- GOCCdown[]
dfgoccdown_simplify <- GOCCdown_simplify[]

plot_goccup <- enrichplot::dotplot(GOCCup, showCategory = 8, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'GOCC up') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_goccdown <- enrichplot::dotplot(GOCCdown, showCategory = 8, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'GOCC down') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_gocc <- ggpubr::ggarrange(plot_goccup, plot_goccdown, nrow = 1, ncol = 2)
ggsave('F5B_DEPs_cancerous_328_GOCC.pdf', plot_gocc, width = 8.5, height = 3.5)

pdf('F5B_DEPs_cancer_up_132_GOCC_similar.pdf', width = 10, height = 8)
similar_matrix_up <- simplifyEnrichment::GO_similarity(dfgoccup$ID, ont = 'CC', db = 'org.Hs.eg.db')
GOCCup_sim <- simplifyEnrichment::simplifyGO(similar_matrix_up)
graphics.off()

pdf('F5B_DEPs_cancer_down_195_GOCC_similar.pdf', width = 10, height = 8)
similar_matrix_down <- simplifyEnrichment::GO_similarity(dfgoccdown$ID, ont = 'CC', db = 'org.Hs.eg.db')
GOCCdown_sim <- simplifyEnrichment::simplifyGO(similar_matrix_down)
graphics.off()

### KEGG ------------
library(KEGG.db)

up_gid <- clusterProfiler::bitr(up_gene, OrgDb = org.Hs.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID')$ENTREZID
down_gid <- clusterProfiler::bitr(down_gene, OrgDb = org.Hs.eg.db, fromType = 'SYMBOL', toType = 'ENTREZID')$ENTREZID

KEGGup <- enrichKEGG(up_gid, organism = 'hsa', use_internal_data = T,
                   keyType = 'ncbi-geneid',
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dfkeggup <- KEGGup[]
KEGGdown <- enrichKEGG(down_gid, organism = 'hsa', use_internal_data = T,
                       keyType = 'ncbi-geneid',
                       minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dfkeggdown <- KEGGdown[]

View(dfkeggup %>% filter(category != 'Human Diseases' | str_detect(Description, '[cC]ancer|[cC]arcino')))
View(dfkeggdown %>% filter(category != 'Human Diseases' | str_detect(Description, '[cC]ancer|[cC]arcino')))

plot_keggup <- KEGGup %>% 
  filter(category != 'Human Diseases' | str_detect(Description, '[cC]ancer|[cC]arcino')) %>% 
  enrichplot::dotplot(orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'KEGG up') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_keggdown <- KEGGdown %>% 
  filter(category != 'Human Diseases' | str_detect(Description, '[cC]ancer|[cC]arcino')) %>% 
  enrichplot::dotplot(orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'KEGG down') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_kegg <- ggpubr::ggarrange(plot_keggup, plot_keggdown, nrow = 1, ncol = 2)
ggsave('F5B_DEPs_cancerous_328_KEGG.pdf', plot_kegg, width = 8.5, height = 3.5)



# # check hsa05210 (Colorectal cancer - Homo sapiens (human))
# library(limma)
# tab <- getGeneKEGGLinks(species = "hsa")
# tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,
#                      column = "SYMBOL", keytype = "ENTREZID")
# tab_ <- tab %>% full_join(DEA, by = c(Symbol = 'Genes'), relationship = 'many-to-many')
# tab_ %>% dplyr::filter(PathwayID == 'hsa05210',
#                        Protein %in% prot_cancer2,
#                        Comparison %in% my_comparisons) %>% View()


### Reactome ------
library(ReactomePA)
# gsePathway()

Reactomeup <- enrichPathway(up_gid,
                            organism = 'human', readable = T,
                            minGSSize = 10, maxGSSize = 2500, pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dfreactomeup <- Reactomeup[]
Reactomedown <- enrichPathway(down_gid,
                              organism = 'human', readable = T,
                              minGSSize = 10, maxGSSize = 2500, pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.05)
dfreactomedown <- Reactomedown[]

plot_reactomeup <- enrichplot::dotplot(Reactomeup, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'Reactome up') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_reactomedown <- enrichplot::dotplot(Reactomedown, orderBy = 'pvalue', x = '-log10(pvalue)', color = 'p.adj', decreasing = F, title = 'Reactome down') +
  labs(x = '-Log10(p value)', y = '') +
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 10, color = 'black'))
plot_reactome <- ggpubr::ggarrange(plot_reactomeup, plot_reactomedown, nrow = 1, ncol = 2)
ggsave('F5B_DEPs_cancerous_328_Reactome.pdf', plot_reactome, width = 10, height = 6)

# view interested pathways
# gList <- DEP_cancer %>%
#   filter(abs(Log2FC) > log2(2)) %>% 
#   left_join(tmp, c(Genes = 'SYMBOL'), relationship = 'many-to-many') %>% 
#   distinct(ENTREZID, Log2FC) %>%
#   drop_na() %>%
#   arrange(Log2FC) %>%
#   distinct(ENTREZID, .keep_all = T)
# 
# viewPathway("Extracellular matrix organization",
#             readable = TRUE,
#             foldChange = gList$Log2FC %>% setNames(gList$ENTREZID))

pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

## 2.6 GSVA ------
# prepare gene matrix
dfgg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_demo_result/N20230320ProteomEx_v2_demo_whole_library_Yueliang/Whole_lib_report.gg_matrix.tsv')
matgg <- dfgg %>% column_to_rownames('Genes') %>% log2() # log2 transform
colnames(matgg) %<>% str_extract('b\\d+_(\\d+|pool)')
dfgg <- matgg %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .) %>% 
  filter(Batch %in% dfpg$Batch)
matgg <- matgg[, dfgg$Batch[!is.na(dfgg$Region)]]
dfgg <- matgg %>% t() %>% as.data.frame() %>% rownames_to_column('Batch') %>% 
  inner_join(info, .)

# Missing ratio cut-off: <20% in at least 3 Regions
dfgg %>% count(Region)

dfgg_miss <- dfgg %>% group_by(Region) %>% select(-(SampleName:BatchHead)) %>%
  summarise_all(function(y) sum(is.na(y)) / length(y)) %>% 
  mutate(Region = as.character(Region))
dfgg_miss %<>% rbind(c(Region = 'ALL', apply(matgg, 1, function(x) sum(is.na(x)) / length(x))))

tmp <- dfgg_miss %>%
  filter(Region != 'ALL') %>% 
  column_to_rownames('Region') %>%
  apply(2, function(y) sum(y < 0.2) >= 3 )
selected_gg <- names(tmp[tmp])
dfgg_miss %<>% select(Region, all_of(selected_gg))
matgg <- matgg[colnames(dfgg_miss)[-1], ]

### missing ratio control -----
mrgg <- apply(matgg, 1, function(x) sum(is.na(x)) / length(x))
plot(density(mrgg))

x <- seq(0.01, 0.99, 0.01)
y <- sapply(x, function(mrgg_cutoff){
  X <- matgg[mrgg < mrgg_cutoff, ]
  return(sum(is.na(X)) / nrow(X) / ncol(X))
})
data.frame(x, y)
plot(x * 100, y * 100, xlab = '% NA ratio cutoff on protein level', ylab = '% whole matggrix NA ratio')
abline(h = 6.5, v = 50)

# dfgg_miss %>% select(Region, all_of(colnames(dfgg_miss)[-1][dfgg_miss[7, -1] < 0.5])) %>% identical(dfgg_miss) # TRUE

matgg <- matgg[colnames(dfgg_miss)[-1], ]
dfgg %<>% select(SampleName:BatchHead, all_of(rownames(matgg)))
dim(matgg) # 4101 123


library(GSVA)
library(GSEABase)
library(enrichplot)
library(limma)
library(pheatmap)
# H: hallmark gene sets (50 gene sets)
# C2_CP: Canonical pathways (3795 gene sets)
# C5_GO_BP: subset of Gene Ontology gene sets (7647 gene sets)
# C6: oncogenic signature gene sets (189 gene sets)
# C7: immunologic signature gene sets (5219 gene sets)
gsH <- getGmt("data/gmt/h.all.v2023.2.Hs.symbols.gmt")
gsC2CP <- getGmt("data/gmt/c2.cp.v2023.2.Hs.symbols.gmt")
gsC5GOBP <- getGmt("data/gmt/c5.go.bp.v2023.2.Hs.symbols.gmt")
gsC6 <- getGmt("data/gmt/c6.all.v2023.2.Hs.symbols.gmt")
gsC7 <- getGmt("data/gmt/c7.all.v2023.2.Hs.symbols.gmt")

### 2.5.1 C5_GO_BP: subset of Gene Ontology gene sets -----



### 2.5.2 H: hallmark gene sets -----
GSVA_h <- GSVA::gsva(
  expr = as.matrix(matgg),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = gsH,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)
identical(colnames(GSVA_h), dfgg$Batch) # TRUE

Type <- as.factor(dfgg$Region)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(
  contrasts = c('TypeL-TypeN', 'TypeH-TypeN', 'TypeC-TypeN', 'TypePC-TypeN', 'TypeCC-TypeN'),
  levels = design
)
fit1 <- lmFit(GSVA_h, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
# res_h <- topTable(fit3, number = Inf)
# qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
# abline(0,1)

res_h <- list(
  L_N = topTable(fit3, coef = 'TypeL-TypeN', number = Inf) %>% rownames_to_column('ID'),
  H_N = topTable(fit3, coef = 'TypeH-TypeN', number = Inf) %>% rownames_to_column('ID'),
  C_N = topTable(fit3, coef = 'TypeC-TypeN', number = Inf) %>% rownames_to_column('ID'),
  PC_N = topTable(fit3, coef = 'TypePC-TypeN', number = Inf) %>% rownames_to_column('ID'),
  CC_N = topTable(fit3, coef = 'TypeCC-TypeN', number = Inf) %>% rownames_to_column('ID')
) %>%
  plyr::ldply(.id = 'Coef')

res_h <- topTable(fit3, number = Inf) %>%
  dplyr::select(AveExpr:adj.P.Val) %>% 
  set_colnames(str_c('Total_', colnames(.))) %>% 
  rownames_to_column('ID') %>% 
  inner_join(res_h, .)

depw_h <- res_h %>% 
  dplyr::filter(Total_adj.P.Val < 0.05) %>%
  plyr::daply('ID', function(dfsub){
    ret <- dfsub %>% filter(adj.P.Val < 0.05)
    ifelse((all(ret$logFC > 0) | all(ret$logFC < 0)) &
             # any(ret$logFC < quantile(res_h$logFC, 0.01) |
             #       ret$logFC > quantile(res_h$logFC, 0.99)) &
             nrow(ret) > 1,
           T, F)
  }) %>% which() %>% names()
res_h %>% filter(ID %in%depw_h) %>% dim()
res_h %>% filter(ID %in%depw_h) %>% count(ID) %>% dim()
x <- res_h %>% filter(ID %in% depw_h) %>% pull(ID) %>% unique()
geneIds(gsH)[x] %>% sapply(length) %>% quantile()
str_replace_all(x, '_', ' ') %>% writeClipboard()
# str_subset(unique(res_h$ID), 'EXTRACELLULAR') %>% .[. %in% x]
# res_h[res_h$ID == 'GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]
# GSVA_h['GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY', ]

res_h %>% filter(ID %in% x) %>% arrange(ID)
X2 <- GSVA_h[x, dfheat$Batch]
rownames(X2) %<>% str_replace_all('_', ' ')
clust_col <- hclust(dist(t(X2), method = 'euclidean'), method = 'ward.D2') %>% 
  as.dendrogram() %>% 
  reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
  as.hclust()
heatmap_GSVA_H <- pheatmap(X2, scale = 'row',
                           annotation_col = dfheat[, 'Region', drop = F],
                           annotation_colors = anno_colors,
                           cluster_rows = T, cluster_cols = F,
                           show_rownames = T, show_colnames = F,
                           color = heat_colors, breaks = bk,
                           fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                           border_color = F, na_col = '#AAAAAA',
                           gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
                           filename = 'F5B_GSVA_Hallmark_heatmap.pdf',
                           width = 7, height = 3.5
)




res_h %>% filter(ID %in%depw_h) %>% arrange(Total_adj.P.Val) %>% pull(ID) %>% unique() %>% head(10) %>% str_replace_all('_', ' ') %>% writeClipboard()







## 2.7 GSEA -------
# prepare gene list
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)

# maxabs <- function(x) x[which.max(abs(x))]
# DEP_cancer2_gsea <- DEP_cancer2 %>% group_by(Protein, Genes) %>% 
#   summarise(Log2FCMax = maxabs(Log2FC),
#             adj.P.min = min(adj.P)) %>%
#   ungroup() %>% 
#   arrange(desc(Log2FCMax))
# tmp1 <- clusterProfiler::bitr(DEP_cancer2_gsea$Genes, OrgDb = 'org.Hs.eg.db', fromType = 'SYMBOL', toType = 'ENTREZID')
# tmp2 <- clusterProfiler::bitr(DEP_cancer2_gsea$Protein, OrgDb = 'org.Hs.eg.db', fromType = 'UNIPROT', toType = 'ENTREZID')
# tmp <- tmp1 %>% full_join(tmp2) %>% setNames(c('Genes', 'ENTREZID', 'Protein'))
# DEP_cancer2_gsea %<>% left_join(tmp)
# 
# gene_list <- drop_na(DEP_cancer2_gsea, Genes)$Log2FCMax %>% setNames(drop_na(DEP_cancer2_gsea, Genes)$Genes)
# gid_list <- drop_na(DEP_cancer2_gsea, ENTREZID)$Log2FCMax %>% setNames(drop_na(DEP_cancer2_gsea, ENTREZID)$ENTREZID)


# use all detected genes (filtered by adj.P <0.05)
DEA_gsea <- DEA %>% 
  filter(Comparison %in% my_comparisons) %>%
  drop_na(adj.P, Log2FC) %>%
  arrange(adj.P) %>%
  plyr::ddply('Protein', function(dfsub){
    ret <- dfsub %>% 
      slice(1) %>% select(Protein, Genes, Log2FC, adj.P)
    ret$Consistent <- all(dfsub$Log2FC > 0) | all(dfsub$Log2FC < 0)
    return(ret)
  })
DEA_gsea %<>% arrange(desc(Log2FC))
DEA_gsea_sig <- DEA_gsea %>% filter(adj.P < 0.05, Consistent)

tmp1 <- clusterProfiler::bitr(DEA_gsea_sig$Genes, OrgDb = 'org.Hs.eg.db', fromType = 'SYMBOL', toType = 'ENTREZID')
tmp2 <- clusterProfiler::bitr(DEA_gsea_sig$Protein, OrgDb = 'org.Hs.eg.db', fromType = 'UNIPROT', toType = 'ENTREZID')
tmp <- tmp1 %>% full_join(tmp2) %>% setNames(c('Genes', 'ENTREZID', 'Protein'))
DEA_gsea_sig %<>% left_join(tmp)

gene_list <- drop_na(DEA_gsea_sig, Genes)$Log2FC %>% setNames(drop_na(DEA_gsea_sig, Genes)$Genes)
gid_list <- drop_na(DEA_gsea_sig, ENTREZID)$Log2FC %>% setNames(drop_na(DEA_gsea_sig, ENTREZID)$ENTREZID)


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

### GOBP -----
gseGOBP <- gseGO(gene_list, OrgDb = org.Hs.eg.db,
                 ont = 'BP', keyType = 'SYMBOL',
                 minGSSize = 10, maxGSSize = 3000,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
df_gsegobp <- gseGOBP[]

# goplot(gseGOBP, showCategory = 10)
# gseaplot2(gseGOBP, 'GO:0090304', pvalue_table = T, ES_geom = 'line')
plot_gseGOBP <- gseGOBP %>%
  filter(p.adjust < 1e-3, abs(NES) > 1.8, setSize < 300) %>%
  ridgeplot(showCategory = 20, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5B_DEPs_cancerous_328_GSEA_GOBP.pdf', plot_gseGOBP, width = 9, height = 9)


### GOCC -----
gseGOCC <- gseGO(gene_list, OrgDb = org.Hs.eg.db,
                 ont = 'CC', keyType = 'SYMBOL',
                 minGSSize = 30, maxGSSize = 5000,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
df_gsegocc <- gseGOCC[]

# goplot(gseGOCC, showCategory = 10)
# gseaplot2(gseGOCC, 'GO:0090304', pvalue_table = T, ES_geom = 'line')
plot_gseGOCC <- gseGOCC %>%
  filter(p.adjust < 1e-3, abs(NES) > 1.8) %>%
  ridgeplot(showCategory = 20, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5B_DEPs_cancerous_328_GSEA_GOCC.pdf', plot_gseGOCC, width = 9, height = 9)


### GOMF -----
gseGOMF <- gseGO(gene_list, OrgDb = org.Hs.eg.db,
                 ont = 'MF', keyType = 'SYMBOL',
                 minGSSize = 30, maxGSSize = 5000,
                 exponent = 1, eps = 0,
                 pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                 verbose = T, seed = F, by = 'fgsea'
)
df_gsegomf <- gseGOMF[]

# goplot(gseGOMF, showCategory = 10)
# gseaplot2(gseGOMF, 'GO:0090304', pvalue_table = T, ES_geom = 'line')
plot_gseGOMF <- gseGOMF %>%
  filter(p.adjust < 1e-3, abs(NES) > 1.8, setSize < 1000) %>%
  ridgeplot(showCategory = 20, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5B_DEPs_cancerous_328_GSEA_GOMF.pdf', plot_gseGOMF, width = 9, height = 9)


### KEGG -----
library(KEGG.db)

gseKEGG <- gseKEGG(gid_list, organism = 'hsa', use_internal_data = T,
                   keyType = 'ncbi-geneid',
                   minGSSize = 5, maxGSSize = 2500,
                   exponent = 1, eps = 0,
                   pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                   verbose = T, seed = F, by = 'fgsea')
df_gsekegg <- gseKEGG[]

# gseaplot2(gseKEGG, 'hsa00190', pvalue_table = T, ES_geom = 'line')
# gseaplot2(gseKEGG, 'hsa04512', pvalue_table = T, ES_geom = 'line')
plot_gseKEGG <- gseKEGG %>%
  filter(!str_detect(Description, 'disease')) %>% 
  filter(p.adjust < 1e-3, abs(NES) > 1.8, setSize < 100) %>%
  ridgeplot(showCategory = 20, orderBy = 'pvalue', decreasing = T)+
  scale_fill_continuous(type = 'gradient') +
  theme_classic() +
  theme(text = element_text(size = 15, color = 'black'))
ggsave('F5B_DEPs_cancerous_328_GSEA_KEGG.pdf', plot_gseKEGG, width = 9, height = 9)



### interesting GSEA pathways ------
gseaplot2(gseGOBP, 'GO:0009060') # aerobic respiration

gseaplot2(gseGOMF, 'GO:0003823') # antigen binding
gseaplot2(gseGOMF, 'GO:0005201') # extracellular matrix structural constituent

# gseaplot2(gseKEGG, 'hsa04115') # p53;
gseaplot2(gseKEGG, 'hsa04110') # cell cycle;
gseaplot2(gseKEGG, 'hsa03030') # DNA replication


plot_gsea_hsa00190 <- enrichplot::gseaplot2(gseKEGG, 'hsa00190', title = 'KEGG: Oxidative phosphorylation', pvalue_table = T, ES_geom = 'line')
plot_gsea_hsa04110 <- enrichplot::gseaplot2(gseKEGG, 'hsa04110', title = 'KEGG: Cell cycle', pvalue_table = T, ES_geom = 'line')
plot_gsea_hsa03030 <- enrichplot::gseaplot2(gseKEGG, 'hsa03030', title = 'KEGG: DNA replication', pvalue_table = T, ES_geom = 'line')

plot_gsea_GO0003823 <- enrichplot::gseaplot2(gseGOMF, 'GO:0003823', title = 'GOMF: Antigen binding', pvalue_table = T, ES_geom = 'line')
plot_gsea_GO0005201 <- enrichplot::gseaplot2(gseGOMF, 'GO:0005201', title = 'GOMF: Extracellular matrix structural constituent', pvalue_table = T, ES_geom = 'line')
plot_gsea_GO0006119 <- enrichplot::gseaplot2(gseGOBP, 'GO:0006119', title = 'GOBP: Oxidative phosphorylation', pvalue_table = T, ES_geom = 'line')
plot_gsea_GO0009060 <- enrichplot::gseaplot2(gseGOBP, 'GO:0009060', title = 'GOBP: Aerobic respiration', pvalue_table = T, ES_geom = 'line')
# plot_gsea_GO0005576 <- enrichplot::gseaplot2(gseGOCC, 'GO:0005576', title = 'GOCC: Extracellular region', pvalue_table = T, ES_geom = 'line')
plot_gsea_GO0070469 <- enrichplot::gseaplot2(gseGOCC, 'GO:0070469', title = 'GOCC: Respirasome', pvalue_table = T, ES_geom = 'line')
plot_gsea_GO0098803 <- enrichplot::gseaplot2(gseGOCC, 'GO:0098803', title = 'GOCC: Aerobic respiration', pvalue_table = T, ES_geom = 'line')





pdf('F5B_cancerous_GSEA_gseaplot.pdf', width = 6, height = 6)
print(plot_gsea_hsa00190)
print(plot_gsea_hsa04110)
print(plot_gsea_hsa03030)
print(plot_gsea_GO0003823)
print(plot_gsea_GO0005201)
print(plot_gsea_GO0006119)
print(plot_gsea_GO0009060)
# print(plot_gsea_GO0005576)
print(plot_gsea_GO0070469)
print(plot_gsea_GO0098803)
graphics.off()



# ### boxplot ------
# library(limma)
# tab <- getGeneKEGGLinks(species = "hsa")
# tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,
#                      column = "SYMBOL", keytype = "ENTREZID")
# tab_ <- tab %>% full_join(DEA, by = c(Symbol = 'Genes'), relationship = 'many-to-many')
# 
# pacman::p_unload(pacman::p_loaded(), character.only = T)
# library(magrittr)
# library(tidyverse)
# library(RColorBrewer)
# 
# 
# #### hsa00190 Oxidative phosphorylation-------
# prot_hsa00190 <- tab_ %>%
#   filter(PathwayID == 'hsa00190', Comparison %in% my_comparisons,
#          adj.P < 0.05, Log2FC < 0) %>% pull(Protein) %>% unique() # 75 proteins
# 
# dfbox <- dfpg %>%
#   filter(!is.na(Region)) %>% 
#   select(SampleName:BatchHead, all_of(prot_hsa00190)) %>% 
#   pivot_longer(cols = -c(SampleName:BatchHead), names_to = 'Protein', values_to = 'Log2Intensity') %>% 
#   left_join(DEA %>% distinct(Protein, Genes, Label), relationship = 'many-to-many') %>%
#   mutate(Protein = factor(Protein, levels = prot_hsa00190))
# 
# boxplot_ls <- plyr::dlply(dfbox, 'Protein', function(dfsub){
#   label <- dfsub$Label[1]
#   set.seed(10)
#   ggplot(dfsub, aes(x = Region, y = `Log2Intensity`, color = Region)) +
#     stat_boxplot(geom = 'errorbar', width = 0.4)+
#     geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
#     geom_jitter(alpha = 0.4, size = 1, width = 0.15)+
#     labs(x = '', y = '', subtitle = label) +
#     scale_color_manual(values = mycolors) +
#     theme_classic() +
#     theme(text = element_text(size = 15), legend.position = 'none') +
#     ggpubr::stat_compare_means(
#       method = 't.test',
#       map_signif_level = F, hide.ns = F, label = 'p.signif',
#       comparisons = str_split(my_comparisons, ' vs. '),
#       size = 5, hjust = 0.5, vjust = 0)
# })
# boxplots <- ggpubr::ggarrange(plotlist = boxplot_ls, nrow = 15, ncol = 5)
# ggsave('F5B_cancerous_GSEA_hsa00190_oxidative_phosphorylation_box.pdf', boxplots, width = 5 * 4, height = 15 * 4, limitsize = F)
# 
# 
# gid_hsa05202 <- tab_ %>%
#   filter(PathwayID == 'hsa05202', Comparison %in% my_comparisons,
#          adj.P < 0.05, Log2FC > 0) %>% pull(GeneID) %>% unique()
# gid_hsa00190 <- tab_ %>%
#   filter(PathwayID == 'hsa00190', Comparison %in% my_comparisons,
#          adj.P < 0.05, Log2FC < 0) %>% pull(GeneID) %>% unique()
# 
# df_KEGGprot <- list(`Transcriptional misregulation in cancer (hsa05202)` = gid_hsa05202,
#                     `Oxidative phosphorylation (hsa00190)` = gid_hsa00190) %>% 
#   plyr::ldply(.id = 'KEGG pathway', .fun = function(gids){
#     fanyi::gene_summary(gids)
#   })
# df_KEGGprot %<>%
#   rename(GeneID = uid) %>%
#   left_join(tab_ %>% filter(Comparison %in% my_comparisons) %>% select(GeneID, Protein, Comparison, Log2FC, P, adj.P, MissingRatio_N, MissingRatio_L, MissingRatio_H, MissingRatio_C, MissingRatio_PC, MissingRatio_CC))
# rio::export(df_KEGGprot, 'F5B_cancerous_GSEA_KEGG_pathway_proteins.xlsx')
# 
# 
# 
# 
# 
# #function to convert entrezids to gene name
# # gids2gnms <- function(entrezid_list) {
# #   # every element contains gids seperated by /
# #   gid_list <- str_split(entrezid_list, '/')
# #   sapply(gid_list, function(gids){
# #     gnms <- clusterProfiler::bitr(gids, 
# #                           OrgDb = 'org.Hs.eg.db', 
# #                           fromType = "ENTREZID", 
# #                           toType = "SYMBOL")$SYMBOL %>% 
# #       str_c(collapse = '/')
# #   })
# # }
# # gseKEGG_tmp <- gseKEGG
# # gseKEGG_tmp@result$core_enrichment <- sapply(gseKEGG_tmp@result$core_enrichment, gids2gnms)
# # gseKEGG_tmp %>% filter(ID == 'hsa00190') %>% cnetplot(foldChange = )
# gseKEGG_tmp <- DOSE::setReadable(gseKEGG, 'org.Hs.eg.db', 'ENTREZID')
# gseKEGG_tmp %>% filter(ID == 'hsa00190') %>% cnetplot(foldChange = gid_alllist)
# 
# 
# 
# ### heatmap ------
# library(pheatmap)
# 
# pm_has00190 <- dfpg_dea %>%
#   arrange(Region, Patient) %>% 
#   select(Batch, all_of(prot_hsa00190)) %>% 
#   # group_by(Region) %>% summarise_all(function(x) log2(mean(2^x, na.rm = T))) %>% 
#   column_to_rownames('Batch') %>% 
#   scale() %>% t() # Z-score
# 
# X <- pm_has00190[intersect(prot_hsa00190, prot_cancer2), ] # 20
# X_labels <- protinfo %>% set_rownames(.$Protein.Ids) %>% .[rownames(X), 'Genes'] %>% str_c(., ' (', rownames(X), ')')
# rownames(X) <- X_labels
# 
# clust_col <- hclust(dist(t(pm_has00190), method = 'euclidean'), method = 'ward.D2') %>%
#   as.dendrogram() %>%
#   reorder(wts = 1:nrow(dfheat)) %>% # reorder to keep original order as much as possible
#   as.hclust()
# heat00190 <- pheatmap(X,
#                       annotation_col = dfheat %>% select(Region, Patient),
#                       annotation_colors = anno_colors,
#                       cluster_rows = T, cluster_cols = clust_col,
#                       show_rownames = T, show_colnames = F,
#                       color = heat_colors, breaks = bk,
#                       fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
#                       border_color = F, na_col = '#AAAAAA',
#                       # gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
#                       filename = 'F5B_DEPs_cancerous_oxidative_phosphorylation_heatmap_cluster.pdf',
#                       width = 6, height = 4.5
# )
# 
# heat00190 <- pheatmap(X,
#                       annotation_col = dfheat %>% select(Region, Patient),
#                       annotation_colors = anno_colors,
#                       cluster_rows = T, cluster_cols = F,
#                       show_rownames = T, show_colnames = F,
#                       color = heat_colors, breaks = bk,
#                       fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
#                       border_color = F, na_col = '#AAAAAA',
#                       gaps_col = cumsum(table(dfheat$Region))[-length(cumsum(table(dfheat$Region)))],
#                       filename = 'F5B_DEPs_cancerous_oxidative_phosphorylation_heatmap.pdf',
#                       width = 6, height = 3.5
# )
# 
# 
# 
# 
# 
# ### Reactome -----
# 
# gseKEGG@result %>% filter(ID == 'hsa00190') %>% pull(core_enrichment) %>% str_split('/')
# 
# 
# 
# 
# 
# gseKEGG@result %>% filter(ID == 'hsa00190') %>% pull(core_enrichment) %>% str_split('/') %>% .[[1]] %>% intersect(names(gid_alllist))
# 
# tab_ %>%
#   filter(PathwayID == 'hsa00190', Comparison %in% my_comparisons,
#          adj.P < 0.05, Log2FC < 0) %>% pull(GeneID) %>% unique() %>% 
#   intersect(names(gid_alllist))
# 



## 2.x output data --------------------
list(DEP_cancer328 = DEP_cancer2,
     DEP_cancer328_expr = dfbox,
     DEP_cancer328_heatmap = X %>% as.data.frame() %>% rownames_to_column(),
     DEP_cancer328_UMAP = df_umap) %>%
  rio::export('F5B_DEP_cancerous.xlsx')

pkgs <- data.frame(Package = c('clusterProfiler', 'ReactomePA', 'org.Hs.eg.db'))
pkgs$Version <- sapply(pkgs$Package, function(pkg) as.character(packageVersion(pkg)))
list(GOBPup_sim_sim = GOBPup_simplify_sim %>% arrange(cluster) %>% inner_join(dfgobpup_simplify, ., by = c(ID = 'id')),
     GOBPdown_sim_sim = GOBPdown_simplify_sim %>% arrange(cluster) %>% inner_join(dfgobpdown_simplify, ., by = c(ID = 'id')),
     GOBPup_similar_mat = similar_matrix_up %>% as.data.frame() %>% rownames_to_column(),
     GOBPdown_similar_mat = similar_matrix_down %>% as.data.frame() %>% rownames_to_column(),
     GOCCup_sim = GOCCup_sim %>% arrange(cluster) %>% inner_join(dfgoccup, ., by = c(ID = 'id')),
     GOCCdown_sim = GOCCdown_sim %>% arrange(cluster) %>% inner_join(dfgoccdown, ., by = c(ID = 'id')),
     GOCCup_similar_mat = similar_matrix_up %>% as.data.frame() %>% rownames_to_column(),
     GOCCdown_similar_mat = similar_matrix_down %>% as.data.frame() %>% rownames_to_column(),
     KEGGup = dfkeggup,
     KEGGdown = dfkeggdown,
     Reactomeup = dfreactomeup,
     Reactomedown = dfreactomedown,
     PackageVersion = pkgs) %>% 
  rio::export('F5B_DEPs_cancerous_328_pathway_enrichment.xlsx')

list(GSEA_GOBP = df_gsegobp,
     GSEA_GOCC =  df_gsegocc,
     GSEA_GOMF = df_gsegomf,
     GSEA_KEGG = df_gsekegg) %>% 
  rio::export('F5B_cancerous_GSEA.xlsx')


save.image('Figure5B_20240108.RData')

# 0.Set environment ----------
#rm(list = ls())
#gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)


# 1.Read data and preprocess --------
# df_lib <- rio::import('data/Mouse_liver_PCT_2g_30_fractions_timsTOF_fragpipe_library_from_GHH.tsv')
load('data/Mouse_liver_PCT_2g_30_fractions_timsTOF_fragpipe_library_from_GHH.RData')
df <- rio::import('data/Hybrid_library_free/Hybrid_library_free_report.pg_matrix.tsv')
colnames(df) %<>% str_split('\\\\') %>% sapply(tail, 1)

df %<>% filter(!str_detect(Protein.Ids, 'CON_')) # remove contaminants
df %>% filter(Protein.Group != Protein.Ids) # 106

df_prot <- df %>% select(Protein.Group:First.Protein.Description) # protein information table
df %<>% select(-Protein.Group, -Protein.Names, -First.Protein.Description) # expression table

mat <- df %>% column_to_rownames('Protein.Ids') %>% select(-Genes) %>% t() %>% log2() # protein matrix (log2-transformed)
df_info <- data.frame(FileName = colnames(df)[-(1:2)],
                      Type = c(rep('Cells', 3), rep('Nuclei', 7)))
df_info$label <- str_extract(df_info$FileName, 'DIA_\\d|\\d\\-\\d') %>% str_replace('^DIA_|^\\d\\-', '') %>% str_c(df_info$Type, .) # sample information table

dfpr <- rio::import('data/Hybrid_library_free/Hybrid_library_free_report.pr_matrix.tsv')
colnames(dfpr) %<>% str_split('\\\\') %>% sapply(tail, 1)
dfpr %<>% select(-Protein.Group, -Protein.Names, -First.Protein.Description)
matpep <- dfpr %>%
  select(Stripped.Sequence, matches('\\.raw$|\\.d$')) %>% 
  group_by(Stripped.Sequence) %>% 
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>% 
  t() %>% log2() # peptide matrix (log2-transformed)


# 2. data analysis --------------
## 2.1 peptides/protein groups (D) ---------
# peptide numbers
dfbar_pep <- apply(matpep, 1, function(x) sum(!is.na(x))) %>% setNames(df_info$label) %>%
  data.frame(Peptide = .) %>% rownames_to_column('label') %>% left_join(df_info)

set.seed(10)
p3D_pep <- ggplot(dfbar_pep, aes(x = Type, y = Peptide, color = Type))+
  stat_boxplot(geom = 'errorbar', width = 0.4)+
  geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '# peptides') +
  scale_y_continuous(labels = scales::scientific, limits = c(0, 1.2e4)) +
  scale_color_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('Cells', 'Nuclei')),
    size = 5, hjust = 0.5, vjust = 0)

# protein group numbers
dfbar_pro <- apply(mat, 1, function(x) sum(!is.na(x))) %>% setNames(df_info$label) %>%
  data.frame(Protein = .) %>% rownames_to_column('label') %>% left_join(df_info)

set.seed(10)
p3D_pro <- ggplot(dfbar_pro, aes(x = Type, y = Protein, color = Type))+
  stat_boxplot(geom = 'errorbar', width = 0.4)+
  geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '# protein groups') +
  scale_y_continuous(limits = c(0, 2200)) +
  scale_color_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('Cells', 'Nuclei')),
    size = 5, hjust = 0.5, vjust = 0)
p3D <- ggpubr::ggarrange(plotlist = list(p3D_pep, p3D_pro), nrow = 1)
ggsave('F3D_number.pdf', p3D, width = 5, height = 5)


## 2.2 venn ---------
df_cell <- df %>% select(1:2, all_of(df_info$FileName[df_info$Type == 'Cells'])) %>%
  filter(!if_all(-(Protein.Ids:Genes), is.na)) # remove all NA; 1933 x 3
df_nucleus <- df %>% select(1:2, all_of(df_info$FileName[df_info$Type == 'Nuclei'])) %>%
  filter(!if_all(-(Protein.Ids:Genes), is.na)) # remove all NA; 951 x 3

prot_cell <- unique(unlist(str_split(df_cell$Protein.Ids, ';')))
prot_nucleus <- unique(unlist(str_split(df_nucleus$Protein.Ids, ';')))
nucleus_specific <- setdiff(prot_nucleus, prot_cell)
cell_specific <- setdiff(prot_cell, prot_nucleus)
nucleus_cell_overlapp <- intersect(prot_cell, prot_nucleus)

prot_list <- list(cells = prot_cell, nuclei = prot_nucleus)
p <- VennDiagram::venn.diagram(x = prot_list,
                               resolution = 300,
                               alpha=rep(0.95, length(prot_list)),
                               # fill=allFills[c(1, 4, 5)],
                               fill = 'white',
                               col = c(Cells = 'green3', Nuclei = 'blue3'),
                               main=stringr::str_glue("Proteins ({length(unique(unlist(prot_list)))} in total)"),
                               #sub = rep,
                               main.cex = 4,
                               sub.cex = 3,
                               cex = 4,
                               cex.lab=4,
                               cat.cex=4,
                               imagetype = "tiff",
                               filename = NULL, disable.logging = T
)
pdf('cells_nuclei_venn.pdf', width = 10, height = 10)
grid::grid.newpage(); grid::grid.draw(p)
graphics.off()

venn <- list(Nucleus_specific_proteins258 = data.frame(Protein = nucleus_specific),
             Cell_specific_protein1276 = data.frame(Protein = cell_specific),
             Nucleus_cell_overlap_protein835 = data.frame(Protien = nucleus_cell_overlapp),
             Cells1933 = df_cell %>% select(Protein.Ids, Genes),
             Nuclei951 = df_nucleus %>% select(Protein.Ids, Genes),
             Overlapped_726groups = df %>% filter(Genes %in% intersect(df_cell$Genes, df_nucleus$Genes)) %>% select(Protein.Ids, Genes),
             Nucleus_specific_230groups = df %>% filter(Genes %in% (setdiff(df_nucleus$Genes, df_cell$Genes) %>% setdiff(''))) %>% select(Protein.Ids, Genes))
rio::export(venn, 'cell_nucleus_protein_compare.xlsx')


## 2.3 boxplot ---------
boxplot(t(mat), main = 'Total 2164 protein groups')
boxplot(na.omit(t(mat)), main = 'Overlapped 144 protein groups')

mat_center <- scale(t(mat), center = apply(mat, 1, median, na.rm = T), scale = F) %>% t()
boxplot(t(mat_center), main = 'Total 2164 protein groups')
boxplot(na.omit(t(mat_center)), main = 'Overlapped 144 protein groups')

mat_qu <- preprocessCore::normalize.quantiles(as.matrix(t(mat)), copy = T, keep.names = T) %>% t()
boxplot(t(mat_qu), main = 'Total 2164 protein groups')
boxplot(na.omit(t(mat_qu)), main = 'Overlapped 144 protein groups')

## 2.4 intraclass CV of proteins -------
cv <- function(x, na.rm = T) sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)

df_cv <- 2^mat %>% as.data.frame() %>% rownames_to_column('FileName') %>%
  inner_join(df_info, .) %>% 
  group_by(Type) %>% select(-FileName) %>% summarise_if(is.numeric, cv)

pdf('cv.pdf', width = 5, height = 5)
boxplot(t(df_cv %>% column_to_rownames('Type')))
df_cv %>% pivot_longer(cols = -Type, names_to = 'Protein.Ids', values_to = 'CV') %>% 
  ggplot(aes(x = Type, y = CV, color = Type, fill = Type))+
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(fill = 'white', color = 'black', width = 0.1) +
  labs(x = '', y = '', subtitle = '# CV of protein groups') +
  scale_color_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  scale_fill_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')
graphics.off()

dfqu_cv <- 2^mat_qu %>% as.data.frame() %>% rownames_to_column('FileName') %>%
  inner_join(df_info, .) %>% 
  group_by(Type) %>% select(-FileName) %>% summarise_all(cv)
boxplot(t(dfqu_cv %>% column_to_rownames('Type')))
dfqu_cv %>% pivot_longer(cols = -Type, names_to = 'Protein.Ids', values_to = 'CV') %>% 
  ggplot(aes(x = Type, y = CV, color = Type, fill = Type))+
  geom_violin(alpha = 0.6, color = NA) +
  geom_boxplot(fill = 'white', color = 'black', width = 0.1) +
  labs(x = '', y = '', subtitle = '# CV of protein groups') +
  scale_color_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  scale_fill_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

## 2.5 Pearson's correlation (E) --------
mat_cor <- mat
rownames(mat_cor)

cors <- mat %>% set_rownames(df_info$label) %>% t() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>% 
  round(2)
cors_cell <- mat %>% set_rownames(df_info$label) %>% .[1:3, ] %>% t() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>% 
  round(2)
cors_nucleus <- mat %>% set_rownames(df_info$label) %>% .[-(1:3), ] %>% t() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>% 
  round(2)

library(corrplot)
pdf('F3E_pearson_correlation.pdf', width = 8, height = 8)
corrplot::corrplot(cors, method = 'square', order = 'original',
                   type = 'full', diag = T,
                   cl.cex = 1.4, tl.pos = "lt", tl.col = "black",
                   tl.cex = 1.5, tl.srt = 45,
                   # col = corrplot::COL1('YlGn', 10),
                   # is.corr = F, col.lim = c(0.8, 1),
                   addCoef.col = 'white', number.cex = 0.9,
                   mar = c(1, 1, 1, 1), title = "Pearson' correlation (Protein group level)"
)
graphics.off()

list(Figure3D = df_info %>% inner_join(dfbar_pep) %>%
       inner_join(dfbar_pro) %>%
       rename(`# peptides` = Peptide, `# protein groups` = Protein),
     Figure3E = cors %>% as.data.frame() %>% rownames_to_column()) %>% 
  rio::export('F3_D-E.xlsx')



# 3.Enrichment analysis ---------------------------------------------------
## 3.1 Pathway enrichment ----------
library(clusterProfiler)
library(org.Mm.eg.db)
### GO ------------
GOCC <- enrichGO(nucleus_specific, OrgDb = org.Mm.eg.db,
                 ont = 'CC', keyType = 'UNIPROT', readable = T,
                 minGSSize = 10, maxGSSize = 2500, pAdjustMethod = 'BH',
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2)
GOCC_simplify <- clusterProfiler::simplify(GOCC, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgocc <- GOCC[]
dfgocc_simplify <- GOCC_simplify[]
enrichplot::dotplot(GOCC)

# similar_matrix <- simplifyEnrichment::GO_similarity(dfgocc_simplify$ID, ont = 'CC', db = 'org.Mm.eg.db')
# GOCC_simplify_sim <- simplifyEnrichment::simplifyGO(similar_matrix)

### KEGG -----
library(KEGG.db)

tmp <- clusterProfiler::bitr_kegg(nucleus_specific, fromType = 'uniprot', toType = 'ncbi-proteinid', organism = 'mmu')
KEGG <- enrichKEGG(tmp$`ncbi-proteinid`, organism = 'mmu',
                   keyType = 'ncbi-proteinid',
                   minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dfkegg <- KEGG[]
enrichplot::dotplot(KEGG)

### Reactome ------
library(ReactomePA)
# gsePathway()

tmp <- clusterProfiler::bitr(nucleus_specific, fromType = 'UNIPROT', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db')
Reactome <- enrichPathway(tmp$ENTREZID,
                          organism = 'mouse', readable = T,
                          minGSSize = 10, maxGSSize = 2500, pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2)
enrichplot::dotplot(Reactome)

# 
# viewPathway("E2F mediated regulation of DNA replication", 
#             readable = TRUE, 
#             foldChange = geneList)

### compareCluster -----
# GOCC
prot_group <- list(nucleus_specific = nucleus_specific,
                   cell_specific = cell_specific,
                   cell_nucleus_overlap = nucleus_cell_overlapp)
compareGOCC <- compareCluster(
  prot_group, fun="enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = 'CC', keyType = 'UNIPROT', readable = T,
  minGSSize = 10, maxGSSize = 2500, pAdjustMethod = 'BH',
  pvalueCutoff = 0.05, qvalueCutoff = 0.2
)
enrichplot::dotplot(compareGOCC)

# KEGG
pid_group <- list(nucleus_specific = clusterProfiler::bitr_kegg(nucleus_specific, fromType = 'uniprot', toType = 'ncbi-proteinid', organism = 'mmu')$`ncbi-proteinid`,
                  cell_specific = clusterProfiler::bitr_kegg(cell_specific, fromType = 'uniprot', toType = 'ncbi-proteinid', organism = 'mmu')$`ncbi-proteinid`,
                  cell_nucleus_overlap = clusterProfiler::bitr_kegg(nucleus_cell_overlapp, fromType = 'uniprot', toType = 'ncbi-proteinid', organism = 'mmu')$`ncbi-proteinid`
)
compareKEGG <- compareCluster(
  pid_group, fun="enrichKEGG",
  organism = 'mmu',
  keyType = 'ncbi-proteinid',
  minGSSize = 10, maxGSSize = 500, pAdjustMethod = 'BH',
  pvalueCutoff = 0.05, qvalueCutoff = 0.2
)
enrichplot::dotplot(compareKEGG)

# Reactome
entrez_group <- list(nucleus_specific = clusterProfiler::bitr(nucleus_specific, fromType = 'UNIPROT', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db')$ENTREZID,
                     cell_specific = clusterProfiler::bitr(cell_specific, fromType = 'UNIPROT', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db')$ENTREZID,
                     cell_nucleus_overlap = clusterProfiler::bitr(nucleus_cell_overlapp, fromType = 'UNIPROT', toType = 'ENTREZID', OrgDb = 'org.Mm.eg.db')$ENTREZID
)
compareReactome <- compareCluster(
  entrez_group, fun="enrichPathway",
  organism = 'mouse', readable = T,
  minGSSize = 10, maxGSSize = 2500, pAdjustMethod = "BH",
  pvalueCutoff = 0.05, qvalueCutoff = 0.2
)
enrichplot::dotplot(compareReactome)

### output -----
pkgs <- data.frame(Package = c('clusterProfiler', 'ReactomePA', 'org.Mm.eg.db'))
pkgs$Version <- sapply(pkgs$Package, function(pkg) as.character(packageVersion(pkg)))

list(nucleus_GOCC = dfgocc,
     nucleus_GOCC_sim = dfgocc_simplify,
     nucleus_KEGG = dfkegg,
     nucleus_Reactome = Reactome[],
     GOCC_compare = compareGOCC[],
     KEGG_compare = compareKEGG[],
     Reactome_compare = compareReactome[],
     nucleusProt258 = data.frame(Protein = nucleus_specific),
     cellProt1276 = data.frame(Protein = cell_specific),
     overlapProt835 = data.frame(Protein = nucleus_cell_overlapp),
     PackageVersion = pkgs) %>% 
  rio::export('nucleus_specific258prot_enrichment.xlsx')

pdf('nucleus_specific258prot_enrichment.pdf', width = 6, height = 4)
enrichplot::dotplot(GOCC, orderBy = 'p.adjust', decreasing = F, showCategory = 10, title = 'GOCC - nucleus')
enrichplot::dotplot(GOCC_simplify, orderBy = 'p.adjust', decreasing = F, showCategory = 8, title = 'GOCC (simplified) - nucleus')
enrichplot::dotplot(KEGG, orderBy = 'p.adjust', decreasing = F, title = 'KEGG - nucleus')
enrichplot::dotplot(Reactome, orderBy = 'p.adjust', decreasing = F, showCategory = 8, title = 'Reactome - nucleus')
graphics.off()

pdf('nucleus258_cell1276_overlap835_compare_enrichment.pdf', width = 8, height = 6)
enrichplot::dotplot(compareGOCC, showCategory = 3, title = 'GOCC comparison')
enrichplot::dotplot(compareKEGG, showCategory = 3, title = 'KEGG comparison')
enrichplot::dotplot(compareReactome, showCategory = 3, title = 'Reactome comparison')
graphics.off()


## 3.2 GSVA --------
### GSVA score ------
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)

mat_g <- df %>% count(Genes) %>% filter(n > 1) %>% anti_join(df, ., by = 'Genes') %>% 
  column_to_rownames('Genes') %>% select(-Protein.Ids) %>% t() %>% log2()
# mat_g[is.na(mat_g)] <- min(mat_g, na.rm = T) + log2(0.8)

mat_g_center <- scale(t(mat_g), center = apply(mat_g, 1, median, na.rm = T), scale = F) %>% t()
# mat_g_center[is.na(mat_g_center)] <- min(mat_g_center, na.rm = T)
boxplot(t(mat_g_center), main = 'Total 2164 genes')
# boxplot(t(scale(mat_g_center)), main = 'Total 2164 genes')
boxplot(t(mat_g_center[, Reduce(intersect, list(df_cell$Genes, df_nucleus$Genes, colnames(mat_g_center)))]), main = 'Intersect 720 genes')

mat_g_qu <- preprocessCore::normalize.quantiles(as.matrix(t(mat_g)), copy = T, keep.names = T) %>% t()
mat_g_qu[is.na(mat_g_qu)] <- min(mat_g_qu, na.rm = T) + log2(0.8)
boxplot(t(mat_g_qu), main = 'Total 2164 genes')
# boxplot(t(scale(mat_g_qu)), main = 'Total 2164 genes')
boxplot(t(mat_g_qu[, Reduce(intersect, list(df_cell$Genes, df_nucleus$Genes, colnames(mat_g_qu)))]), main = 'Intersect 720 genes')


library(GSVA)
library(GSEABase)

g5mM <- getGmt("data/gmt/m5.all.v2023.2.Mm.symbols.gmt")
g2mM <- getGmt("data/gmt/m2.all.v2023.2.Mm.symbols.gmt")
GSVA_g5 <- GSVA::gsva(
  expr = as.matrix(t(mat_g_qu)),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = g5mM,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)
GSVA_g2 <- GSVA::gsva(
  expr = as.matrix(t(mat_g_qu)),# a matrix of expression values where rows correspond to genes and columns correspond to samples.
  gset.idx.list = g2mM,
  method = 'gsva', kcdf = 'Gaussian', abs.ranking = T
)

pheatmap::pheatmap(
  GSVA_g5[!apply(GSVA_g5, 1, function(x) all(x == 1)), ],
  # breaks = my_breaks, color = my_colors,
  color = colorRampPalette(c("blue", "white","red"))(100),
  scale = 'none', cluster_rows = T, cluster_cols = T,
  annotation_col = df_info %>% column_to_rownames('FileName'),
  # annotation_row = ann_row,
  # annotation_colors = list(),
  show_rownames = F,
  show_colnames = T, angle_col = 45,
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'M5: ontology gene sets GSVA'
)
pheatmap::pheatmap(
  GSVA_g2[!apply(GSVA_g2, 1, function(x) all(x == 1)), ],
  # breaks = my_breaks, color = my_colors,
  color = colorRampPalette(c("blue", "white","red"))(100),
  scale = 'none', cluster_rows = T, cluster_cols = T,
  annotation_col = df_info %>% column_to_rownames('FileName'),
  # annotation_row = ann_row,
  # annotation_colors = list(),
  show_rownames = F,
  show_colnames = T, angle_col = 45,
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = 'M2: curated gene sets'
)

### DEA ------
library(limma)

identical(colnames(GSVA_g5), df_info$FileName) # TRUE
identical(colnames(GSVA_g2), df_info$FileName) # TRUE

Type <- as.factor(df_info$Type)
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(Type = 'TypeNuclei-TypeCells', levels = design)

fit1 <- lmFit(GSVA_g5, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
res_g5 <- topTable(fit3, coef = 'Type', number = Inf)
# colnames(res_g2) %<>% str_c('_g2')

fit1 <- lmFit(GSVA_g2, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
res_g2 <- topTable(fit3, coef = 'Type', number = Inf)
qqt(fit3$t, df = fit3$df.prior+fit3$df.residual, pch = 16, cex = 0.2)
abline(0,1)

res_g5 %>% filter(adj.P.Val < 0.05) %>% pull(adj.P.Val) %>% hist()
res_g5 %>% filter(adj.P.Val < 0.05) %>% pull(logFC) %>% hist()
res_g5 %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.2)
res_g5 %>% filter(adj.P.Val < 0.01, abs(logFC) > 0.3)
res_g5 %>% filter(adj.P.Val < 0.01, abs(logFC) > 0.2) %>% rownames_to_column('ID') %>% filter(str_detect(ID, '^GOCC'))


res_g2 %>% filter(P.Value < 0.01) %>% pull(adj.P.Val) %>% hist()
res_g2 %>% filter(P.Value < 0.01) %>% pull(logFC) %>% hist()
res_g2 %>% filter(P.Value < 0.01, abs(logFC) > 0.2)
res_g2 %>% filter(P.Value < 0.01, abs(logFC) > 0.3)

list(`mouse_ontoloty_GSVA` = GSVA_g5 %>% as.data.frame() %>% rownames_to_column('ID'),
     `mouse_ontoloty_Nuclei-Cell` = res_g5 %>% rownames_to_column('ID'),
     `mouse_curated_GSVA` = GSVA_g2 %>% as.data.frame() %>% rownames_to_column('ID'),
     `mouse_curated_Nuclei-Cell` = res_g2 %>% rownames_to_column('ID')
) %>% rio::export('nucleus_vs_cell_GSVA.xlsx')

## Subcellular location -------
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)

allprot <- union(prot_cell, prot_nucleus)
'' %in% allprot # FALSE
str_subset(allprot, '^CON') # character(0)
writeClipboard(allprot)

df_loc <- rio::import('data/idmapping_2023_11_18.xlsx')
# loc_ls <- lapply(df_loc$`Subcellular location [CC]`, function(x){
#   x %>% str_split('SUBCELLULAR LOCATION: (\\[.+?\\]: )?') %>% .[[1]] %>% 
#     str_split('\\. ') %>% unlist() %>% str_extract('^.+?\\{') %>% str_remove(' \\{$') %>%
#     unique() %>% .[!is.na(.)] %>% .[. != '']
# })
loc_ls <- lapply(df_loc$`Subcellular location [CC]`, function(x){
  x %>% str_split('SUBCELLULAR LOCATION: (\\[.+?\\]: )?') %>% .[[1]] %>% 
    str_split('\\. ') %>% unlist() %>% str_extract('^.+?\\{|^[\\w, ]+$') %>% str_remove(' \\{$') %>%
    unique() %>% .[!is.na(.)] %>% .[. != '']
})
loc_ls_simple <- lapply(loc_ls, function(x){
  str_extract(x, '^[\\w ]+') %>% unique() %>% sort()
})

nucleus_loc <- c('Nucleus', 'Chromosome', 'Nucleus membrane', 'Nucleus envelope', 'Nucleus speckle', 'Nucleus lamina', 'Nucleus matrix', 'Nucleus inner membrane', 'Nucleus outer membrane',
                 'Cytoplasm')
df_loc1 <- sapply(loc_ls_simple, function(x){ # nucleus-located proteins
  nucleus_loc %in% x
}) %>%
  set_rownames(nucleus_loc) %>% set_colnames(df_loc$From) %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column('From') %>% 
  inner_join(df_loc, .)

df_loc2 <- lapply(seq_along(nucleus_loc), function(i){ # mean of nucleus-located proteins
  loc <- nucleus_loc[i]
  loc_prots <- df_loc1$From[df_loc1[[loc]]]
  prot_ptn <- str_c(loc_prots, collapse = '|')
  mat_tmp <- mat_qu[, str_which(colnames(mat_qu), prot_ptn), drop = F]
  rowMeans(2 ^ mat_tmp, na.rm = T) %>% as.data.frame() %>% setNames(loc)
}) %>% cbind.data.frame() %>% rownames_to_column('FileName') %>% 
  left_join(df_info, .)

set.seed(2023)
df_loc2[is.na(df_loc2)] <- df_loc2 %>% select(Nucleus:`Nucleus outer membrane`) %>% min(na.rm = T) * rnorm(1, 1, 0.01)

df_loc2_long <- df_loc2 %>%
  pivot_longer(cols = -(FileName:label), names_to = 'Location', values_to = 'Intensity') %>%
  mutate(Log2Intensity = log2(Intensity))

df_pvalue <- plyr::ddply(df_loc2_long, 'Location', function(dfsub){
  xn <- dfsub$Log2Intensity[dfsub$Type == 'Nuclei']
  xc <- dfsub$Log2Intensity[dfsub$Type == 'Cells']
  t.test(xn, xc, paired = F, var.equal = F)$p.value %>% setNames('p.value')
})

set.seed(10)
p <- ggplot(df_loc2_long, aes(x = Location, y = Log2Intensity, color = Type))+
  facet_grid(Location~., scale='free') +
  # geom_violin(fill = 'grey', color = NA) +
  # geom_boxplot(width = 1) +
  # stat_boxplot(geom = 'errorbar', width = 0.4)+
  geom_boxplot(fill = '#FFFFFF', width = 0.5, outlier.shape = NA)+
  geom_jitter(alpha = 0.7, size = 3, width = 0.15)+
  coord_flip() +
  labs(x = '') +
  scale_color_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

# ggpubr::geom_signif(comparisons = c('Cell', 'Nucleus'), step_increase = 0.1, map_signif_level = F, test = t.test)
ggsave('subcellular_location.pdf', p, width = 7, height = 4)


list(nucleus_location_protein = df_loc1, nucleus_location_mean = df_loc2,
     boxplot = df_loc2_long, t_test = df_pvalue) %>% 
  rio::export('subcellular_location.xlsx')


# fold-change curve
nafill <- min(mat_qu, na.rm = T)
mat_qu_imp <- mat_qu
set.seed(2023)
mat_qu_imp[is.na(mat_qu_imp)] <- min(mat_qu_imp, na.rm = T) * rnorm(sum(is.na(mat_qu_imp)), 1, 0.01)
df_fc <- lapply(seq_along(nucleus_loc), function(i){ # mean of nucleus-located proteins
  loc <- nucleus_loc[i]
  loc_prots <- df_loc1$From[df_loc1[[loc]]]
  prot_ptn <- str_c(loc_prots, collapse = '|')
  mat_tmp <- mat_qu_imp[, str_which(colnames(mat_qu), prot_ptn), drop = F]
  data.frame(
    loc,
    apply(mat_tmp, 2, function(y) {
      log2(mean(2^y[4:10]) / mean(2^y[1:3])) # nucleus vs. cell
    }),
    apply(mat_tmp, 2, function(y) {
      t.test(y[1:3], y[4:10])$p.value
    })
  ) %>% 
    setNames(c('Location', 'Log2FC', 'p.value')) %>%
    rownames_to_column('Protein.Ids')
}) %>% plyr::ldply() %>% inner_join(df_prot, .)
df_fc$adj.P <- p.adjust(df_fc$p.value, 'BH')
rio::export(df_fc, 'subcellular_location_nucleus-cell_difference.xlsx')

library(ggridges)
pdf('subcellular_location_Log2FC.pdf', width = 7, height = 4)
ggplot(df_fc, aes(x = Log2FC, y = Location, fill = Location)) +#数据
  geom_density_ridges(stat = "density_ridges") +
  labs(x = 'Log2FC(nuclei/cells)', y = '') +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = 'none')

ggplot(df_fc %>% filter(adj.P < 0.05), aes(x = Log2FC, y = Location, fill = Location)) +#数据
  geom_density_ridges(stat = "density_ridges") +
  labs(x = 'Log2FC(nuclei/cells) (adj.P <0.05)', y = '') +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = 'none')
graphics.off()

# nucleus protein ratio
prot_cell <- unique(unlist(str_split(df_cell$Protein.Ids, ';')))
prot_nucleus <- unique(unlist(str_split(df_nucleus$Protein.Ids, ';')))
nucleus_specific <- setdiff(prot_nucleus, prot_cell) # 258
cell_specific <- setdiff(prot_cell, prot_nucleus) # 1276
cn_overlap <- intersect(prot_cell, prot_nucleus) # 835

tmp1 <- df_loc1 %>% filter(Nucleus) %>% pull(From) # Nucleus-located proteins
Nucleus <- c(sum(nucleus_specific %in% tmp1) * 100 / length(nucleus_specific), # 45.7%
             sum(cell_specific %in% tmp1) * 100 / length(cell_specific), # 14.7%
             sum(cn_overlap %in% tmp1) * 100 / length(cn_overlap)) # 36.3%

tmp2 <- df_loc1 %>% filter(Cytoplasm) %>% pull(From) # Nucleus-located proteins
Cytoplasm <- c(sum(nucleus_specific %in% tmp2) * 100 / length(nucleus_specific), # 39.5%
               sum(cell_specific %in% tmp2) * 100 / length(cell_specific), # 32.4%
               sum(cn_overlap %in% tmp2) * 100 / length(cn_overlap)) # 53.8%

df_ratio <- data.frame(Nucleus, Cytoplasm) %>% round(2) %>%
  set_rownames(c('nuclei-specific (%)', 'cells-specific (%)', 'nuclei-cells overlap (%)')) %>%
  t() %>% as.data.frame() %>% 
  rownames_to_column('Subcellular location')
rio::export(df_ratio, 'subcellular_location_nucleus-cell_ratio.xlsx')

p <- df_ratio %>% pivot_longer(-`Subcellular location`) %>% 
  filter(`Subcellular location` == 'Nucleus') %>% 
  ggplot()+
  geom_line(aes(x = name, y = value, group = `Subcellular location`, color = `Subcellular location`), size = 2, color = 'black') +
  geom_label(aes(x = name, y = value, label = value)) +
  labs(x = '', y = '% nuclei') +
  coord_flip() +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = 'none')
ggsave('subcellular_location_nucleus-cell_ratio.pdf', p, width = 6, height = 4)

## DEG analysis ------
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)

mat_g <- df %>% count(Genes) %>% filter(n > 1) %>% anti_join(df, ., by = 'Genes') %>% 
  column_to_rownames('Genes') %>% select(-Protein.Ids) %>% t() %>% log2()
# mat_g[is.na(mat_g)] <- min(mat_g, na.rm = T)
mat_g_qu <- preprocessCore::normalize.quantiles(as.matrix(t(mat_g)), copy = T, keep.names = T) %>% t()
mat_g_qu[is.na(mat_g_qu)] <- min(mat_g_qu, na.rm = T) + log2(0.8)

pheatmap::pheatmap(
  t(mat_g_qu),
  color = colorRampPalette(c("blue", "white","red"))(100),
  scale = 'row', cluster_rows = T, cluster_cols = T,
  annotation_col = df_info %>% column_to_rownames('FileName'),
  # annotation_row = ann_row,
  # annotation_colors = list(),
  show_rownames = F,
  show_colnames = T, angle_col = 45,
  # cellwidth = 10, cellheight = 15,
  fontsize = 10,
  main = ''
)


## UMAP ---------
library(umap)
umap <- umap(mat_g_qu, n_neighbors = 5)
df_umap <- df_info %>% inner_join(umap$layout %>% as.data.frame() %>% rownames_to_column('FileName'))
p_umap <- ggplot(df_umap, aes(x=V1, y=V2,fill=Type,color=Type)) +
  geom_point(size=5,alpha=1)+
  # stat_ellipse(type = "norm",level = 0.95,size=2.5)+
  labs(x="UMAP1",y="UMAP2") +
  scale_color_manual(values = c(Cells = 'green3', Nuclei = 'blue3')) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

save.image('Figure3.RData')

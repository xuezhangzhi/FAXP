# 0.Set environment ----------
#rm(list = ls())
#gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)


# 1.B-J -------------------------------------------------------------------
## 1.1 load data ----
#peptide matrix
dfpep <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx_v2_InTip_vs_InGel_DDA/combined_peptide.tsv')
dfpep %<>% filter(!str_detect(Protein, '^CON_'))

matpep <- dfpep %>%
  select(`Peptide Sequence`, matches('^In\\w+?_1_\\d Intensity$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(matpep) %<>% str_remove(' Intensity$')

matpep_count <- dfpep %>%
  select(`Peptide Sequence`, matches('^In\\w+?_1_\\d.+Spectral Count$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(matpep_count) %<>% str_remove(' Spectral Count$')

#sample info
info <- rio::import('F2_info.xlsx', sheet = 1)
info <- data.frame(SampleName = colnames(matpep),
                   'Label in file' = str_remove(colnames(matpep), '\\d+$'),
                   check.names = F) %>% left_join(info)

#protein matrix
dfpro <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx_v2_InTip_vs_InGel_DDA/combined_protein.tsv')
dfpro %<>% filter(!str_detect(Protein, '^CON_'))
matpro <- dfpro %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^In\\w+?_1_\\d Intensity$')) %>%
  column_to_rownames('Protein ID')
colnames(matpro) %<>% str_remove(' Intensity$')

matpro_count <- dfpro %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^In\\w+?_1_\\d Total Spectral Count$')) %>%
  column_to_rownames('Protein ID')
colnames(matpro_count) %<>% str_remove(' Total Spectral Count$')

#regenerate dfpep and dfpro; within protein groups
dfpep <- matpep %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro <- matpro %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpep_count <- matpep_count %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro_count <- matpro_count %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)

## 1.2 peptide number (B) ---------
dflong <- dfpep_count %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA))
tbl1 <- dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`# Peptides` = sum(!is.na(Count))) %>% 
  left_join(info)

set.seed(0)
p1 <- ggplot(dfbar, aes(x = `Label in figure`, y = `# Peptides`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '# Peptides')+
  scale_color_brewer(palette = 'Set1', direction = -1) +
  scale_y_continuous(labels = scales::scientific, limits = c(0, 2.5e4)) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)

## 1.3 protein number (B) ---------
dflong <- dfpro_count %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA)) %>% 
  separate_rows('Protein')
tbl2 <- dfbar <- dflong %>% 
  group_by(SampleName) %>% 
  summarise(`# Proteins` = sum(!is.na(Count))) %>% 
  left_join(info)

set.seed(0)
p2 <- ggplot(dfbar, aes(x = `Label in figure`, y = `# Proteins`, color = `Label in figure`))+
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '# Proteins')+
  scale_y_continuous(limits = c(0, 3300)) +
  scale_color_brewer(palette = 'Set1', direction = -1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)

## 1.4 Venn (C) --------
pep_list <- plyr::dlply(dfpep_count, '`Label in figure`', function(dfsub){
  dfsub %>% select(-(SampleName:Note)) %>% 
    .[, apply(., 2, function(y) !all(y == 0))] %>% 
    colnames()
})
pro_list <- plyr::dlply(dfpro_count, '`Label in figure`', function(dfsub){
  dfsub %>% select(-(SampleName:Note)) %>% 
    .[, apply(., 2, function(y) !all(y == 0))] %>% 
    colnames() %>% str_split(', ') %>% unlist() %>% unique()
})
# VennDiagram
venn1 <- VennDiagram::venn.diagram(x = pep_list,
                                resolution = 300,
                                alpha=rep(0.95, length(pep_list)),
                                # fill=allFills[c(1, 4, 5)],
                                fill = 'white', col = c(InGel = '#377EB8', InTip = '#E41A1C'),
                                main = stringr::str_glue("Peptides ({nrow(matpep)} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
venn2 <- VennDiagram::venn.diagram(x = pro_list,
                                resolution = 300,
                                alpha=rep(0.95, length(pro_list)),
                                fill = 'white', col = c(InGel = '#E41A1C', InTip = '#377EB8'),
                                main=stringr::str_glue("Proteins ({nrow(matpro)} in total)"),
                                #sub = rep,
                                main.cex = 4,
                                sub.cex = 3,
                                cex = 4,
                                cex.lab=4,
                                cat.cex=4,
                                imagetype = "tiff",
                                filename = NULL, disable.logging = T
)
pdf('F2C_VennDiagram.pdf', width = 10, height = 10)
grid::grid.newpage(); grid::grid.draw(venn1)
grid::grid.newpage(); grid::grid.draw(venn2)
graphics.off()


vennlist <- list(InGel.Peptide = pep_list$InGel,
                 InTip.Peptide = pep_list$InTip,
                 InGel.Protein = pro_list$InGel,
                 InTip.Protein = pro_list$InTip)
venntbl <- matrix(nrow = max(sapply(vennlist, length)), ncol = 4) %>% 
  as.data.frame() %>% setNames(names(vennlist))
for(j in 1:ncol(venntbl)){
  venntbl[1:length(vennlist[[j]]), j] <- vennlist[[j]]
}
rio::export(list(Figure2C = venntbl), 'F2C_VennDiagram.xlsx')

## 1.5 Peptide length curve (D) --------
df <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx_v2_InTip_vs_InGel_DDA/combined_peptide.tsv')
df %<>% filter(!str_detect(`Protein`, '^CON_')) %>% select(`Peptide Sequence`, `Peptide Length`) %>% rename(Peptide = `Peptide Sequence`)

dflong <- dfpep_count %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Count') %>% 
  mutate(Count = ifelse(Count != 0, Count, NA))
tbl3 <- dfbar <- dflong %>% 
  left_join(df) %>% 
  filter(!is.na(Count)) %>% 
  select(-SampleName, -Note, -Count) %>% 
  distinct() %>% 
  count(`Label in figure`, `Peptide Length`)

p3 <- ggplot(dfbar, aes(x = `Peptide Length`, y = n,
                        #color = `Label in figure`,
                        fill = `Label in figure`)) +
  geom_col(position = 'dodge') +
  # geom_point(size = 2) +
  # geom_line() +
  labs(x = '', y = '', subtitle = 'Length of peptides')+
  scale_y_continuous(labels = scales::scientific) +
  scale_color_brewer(palette = 'Set1', direction = -1) +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')


### output B,D-------
pBD <- ggpubr::ggarrange(p1, p2, p3, nrow = 1, widths = c(1, 1, 2))
ggsave('F2_B_D.pdf', pBD, width = 10, height = 5)
tblBD <- list(Figure2B = info %>% inner_join(tbl1) %>% inner_join(tbl2) %>% select(-Note),
     Figure2D = tbl3)
rio::export(tblBD, 'F2_B_D.xlsx')

## 1.6 Physicochemical property (E) --------
library(Peptides)
calc_pep1 <- function(pepseq){
  data.frame(Peptide = pepseq,
             pI_Lehninger = pI(pepseq, pKscale = 'Lehninger'),
             NetCharge_Lehninger7 = pI(pepseq, pKscale = 'EMBOSS'),
             instaIndex = instaIndex(pepseq),
             aIndex = aIndex(pepseq)
  )
}
calc_pep2 <- function(pepseq){
  X1 <- crucianiProperties(pepseq) %>% setNames(pepseq) %>% as.data.frame()
  X2 <- kideraFactors(pepseq) %>% setNames(pepseq) %>% as.data.frame()
  X3 <- zScales(pepseq) %>% setNames(pepseq) %>% as.data.frame()
  rbind(X1, X2, X3) %>% t() %>% as.data.frame() %>% rownames_to_column('Peptide')
}

dfscore1 <- plyr::ldply(pep_list, calc_pep1, .id = 'Label in figure')
dfscore2 <- plyr::ldply(pep_list, calc_pep2, .id = 'Label in figure')
tbl4 <- dfscore <- dfscore1 %>% inner_join(dfscore2)
# maybe too large sample size for KS test
# x1 <- dfscore$Z4[dfscore$`Label in figure` == 'InGel']
# x2 <- dfscore$Z4[dfscore$`Label in figure` == 'InTip']
# ks.test(x1, x2)$p.value <<< 0.05  ----> not in the same distribution

p4_list <- list()
for(i in 3:ncol(dfscore)){
  dftmp <- dfscore %>% select(1, all_of(i)) %>% setNames(c('Label in figure', 'Index'))
  p4_list[[i-2]] <- ggplot(dftmp, aes(x = Index, color = `Label in figure`)) +
    geom_density(alpha = 0.9, linewidth = 1) +
    labs(x = '', y = '', subtitle = str_c('Density of ', colnames(dfscore)[i]))+
    scale_color_brewer(palette = 'Set1', direction = -1) +
    theme_classic() +
    theme(text = element_text(size = 15), legend.position = 'none')
}
p4 <- ggpubr::ggarrange(plotlist = p4_list, nrow = 5, ncol = 5)
ggsave('F2E_identified_peptides_physicochemical_property.pdf', p4, width = 14, height = 14)
rio::export(list(Figure2E = tbl4), 'F2E_identified_peptides_physicochemical_property.xlsx')
# PP1: Polarity,
# PP2: Hydrophobicity,
# PP3: H-bonding
# KF1: Helix/bend preference,
# KF2: Side-chain size,
# KF3: Extended structure preference,
# KF4: Hydrophobicity,
# KF5: Double-bend preference,
# KF6: Partial specific volume,
# KF7: Flat extended preference,
# KF8: Occurrence in alpha region,
# KF9: pK-C,
# KF10: Surrounding hydrophobicity
# Z1: Lipophilicity
# Z2: Steric properties (Steric bulk/Polarizability)
# Z3: Electronic properties (Polarity / Charge)
# Z4 and Z5: They relate electronegativity, heat of formation, electrophilicity and hardness.

## 1.7 Missed cleavage (%) (F) --------
calc_missCleav <- function(pepseq){
  df_missCleav <- data.frame(pepseq = pepseq)
  df_missCleav$MissedCleavage <- 0
  df_missCleav$MissedCleavage <- apply(df_missCleav, 1, function(row){
    Seq <- row['pepseq']
    # exclude the last cleavage site
    if (stringr::str_sub(Seq, -1, -1) == 'P' & stringr::str_sub(Seq, -2, -2) %in% c('K', 'R')){
      Seq <- stringr::str_sub(Seq, 1, -3)
    }else{
      Seq <- stringr::str_sub(Seq, 1, -2)
    }
    #
    Seqs <- stringr::str_sub(Seq, 1:nchar(Seq), 1:nchar(Seq))
    is_KR <- Seqs %in% c('K', 'R')
    isnt_P <- Seqs != 'P'
    flag <- c(isnt_P[2:length(isnt_P)], T)
    mis_cleav <- is_KR & flag
    return(sum(mis_cleav))
  })
  return(df_missCleav)
}

pepList <- dfpep_count %>% 
  select(-Note) %>% 
  plyr::dlply('SampleName', function(dfsub){
  dfsub[, which(dfsub != 0)] %>% select(-(SampleName:'Label in figure')) %>% colnames()
})
df_missCleav <- plyr::ldply(pepList, calc_missCleav, .id = 'SampleName')
tbl5 <- dfbar <- df_missCleav %>%
  count(SampleName, MissedCleavage) %>%
  with_groups(SampleName, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  inner_join(info)

p5 <- ggplot(dfbar, aes(x = `Label in figure`, y = `ratio (%)`, color = `Label in figure`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage))) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  scale_color_brewer(palette = 'Set1', direction = -1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF')) +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)
ggsave('F2F_missed_cleavage.pdf', p5, width = 8, height = 5)

### output B,D,F-------
tblBDF <- list(Figure2B = info %>% inner_join(tbl1) %>% inner_join(tbl2) %>% select(-Note),
              Figure2D = tbl3,
              Figure2F = info %>% inner_join(tbl5) %>% select(-Note))
rio::export(tblBDF, 'F2_BDF.xlsx')


## 1.8 Log2Intensity rank, protein group level (G) --------
#Protein info are downloaded from UniProt.Org, and selected with "liver" detected in "Tissue specificity"
protinfo <- rio::import('uniprotkb_mouse_reviewed_17184_20231213.xlsx') %>% 
  filter(Entry %in% unique(unlist(str_split(colnames(dfpro)[-(1:4)], ', ')))) %>%
  rename(Protein = Entry) %>%
  select(-Reviewed, -Organism)
# mutate(Marker = str_detect(`Tissue specificity`, 'Abundantly expressed in (the )?liver|Enriched in (the )?liver|Expressed abundantly in (the )?liver|Expressed at high levels in (the )?liver|Highest expression in (the )?liver|Highest level in (the )?liver|Liver specific|Liver-specific|Mainly expressed in (the )?liver|Mainly expressed in (the )?liver|Mainly in (the )?liver|Predominantly expressed in (the )?liver|Present at high level in (the )?liver|Strong expression in (the )?liver|Strongly expressed in (the )?liver'))
liverprotinfo <- protinfo %>%
  filter(str_detect(str_to_lower(`Tissue specificity`), 'liver'))

dfrank <- dfpro %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Intensity') %>% 
  mutate(Intensity = ifelse(Intensity != 0, Intensity, NA)) %>% 
  filter(!is.na(Intensity)) %>% 
  group_by(`Label in figure`, Protein) %>% 
  summarise(MeanIntensity = mean(Intensity, na.rm = T)) %>% 
  mutate(Log2MeanIntensity = log2(MeanIntensity)) %>% 
  arrange(desc(Log2MeanIntensity)) %>% 
  with_groups(`Label in figure`, mutate, Rank = 1:length(MeanIntensity))

tbl6 <- dfrank_label <- dfrank %>%
  filter(str_detect(Protein, str_c(liverprotinfo$Protein, collapse = '|')),
         !str_detect(Protein, ', ')) %>%  # do not label protein groups
  left_join(liverprotinfo)

p6 <- ggplot(dfrank, aes(x = Rank, y = Log2MeanIntensity, color = `Label in figure`, fill = `Label in figure`))+
  geom_line(color = '#000000', linewidth = 1) +
  geom_point(data = dfrank_label, size = 1, alpha = 0.5) +
  labs(x = '', y = 'Log2 mean intensity', subtitle = 'Protein group')+
  # scale_y_continuous(limits = c(0, 26)) +
  scale_color_brewer(palette = 'Set1', direction = -1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

p6_new <- ggplot(dfrank, aes(x = Rank, y = Log2MeanIntensity, color = `Label in figure`, fill = `Label in figure`))+
  geom_line(color = '#000000', linewidth = 1) +
  geom_point(data = dfrank_label, size = 1, alpha = 0.5) +
  labs(x = '', y = 'Log2 mean intensity', subtitle = 'Protein group')+
  # scale_y_continuous(limits = c(0, 26)) +
  # scale_color_brewer(palette = 'Set1', direction = -1) +
  scale_color_manual(values = c(InGel = "#00008B", InTip = "#C32022")) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')
ggsave('F2I_intensity_range.pdf', p6_new, width = 5, height = 5)


## 1.9 CV (H) and Pearson's correlation (I) --------
### protein CV -----
dfcv_pro <- dfpro %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Protein', values_to = 'Intensity') %>% 
  mutate(Intensity = ifelse(Intensity != 0, Intensity, NA)) %>% 
  filter(!is.na(Intensity)) %>% 
  group_by(`Label in figure`, Protein) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean)

plot_procv <- ggplot(dfcv_pro, aes(x = `Label in figure`, y = CV, color = `Label in figure`, fill = `Label in figure`))+
  geom_violin(color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Protein group')+
  scale_color_brewer(palette = 'Set1', direction = -1) +
  scale_fill_brewer(palette = 'Set1', direction = -1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)
x1 <- dfcv_pro %>% filter(`Label in figure` == 'InGel') %>% pull(CV)
x2 <- dfcv_pro %>% filter(`Label in figure` == 'InTip') %>% pull(CV)
lbl <- t.test(x1, x2)$p.value # 6.23e-43

### peptide CV -----
dfcv_pep <- dfpep %>%
  pivot_longer(cols = -(SampleName:Note),
               names_to = 'Peptide', values_to = 'Intensity') %>% 
  mutate(Intensity = ifelse(Intensity != 0, Intensity, NA)) %>% 
  filter(!is.na(Intensity)) %>% 
  group_by(`Label in figure`, Peptide) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean)

plot_pepcv <- ggplot(dfcv_pep, aes(x = `Label in figure`, y = CV, color = `Label in figure`, fill = `Label in figure`))+
  geom_violin(color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Peptide')+
  scale_color_brewer(palette = 'Set1', direction = -1) +
  scale_fill_brewer(palette = 'Set1', direction = -1) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(c('InGel', 'InTip')),
    size = 5, hjust = 0.5, vjust = 0)
x1 <- dfcv_pep %>% filter(`Label in figure` == 'InGel') %>% pull(CV)
x2 <- dfcv_pep %>% filter(`Label in figure` == 'InTip') %>% pull(CV)
lbl <- t.test(x1, x2)$p.value # 0

### output G,H -------
pGH <- ggpubr::ggarrange(p6, plot_pepcv, plot_procv, nrow = 1, widths = c(2, 1, 1))
ggsave('F2_G_H.pdf', pGH, width = 10, height = 5)


## 1.10 Pearson's correlation (I) -----
matpro_ <- matpro
matpro_[matpro_ == 0] <- NA
cors_pro <- matpro_ %>% log2() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>%
  round(2)
min(cors_pro) # 0.84

matpep_ <- matpep
matpep_[matpep_ == 0] <- NA
cors_pep <- matpep_ %>% log2() %>%
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>%
  round(2)
min(cors_pep) # 0.56

pdf('F2I_pearson_correlation.pdf', width = 8, height = 8)
corrplot::corrplot(cors_pro, method = 'square', order = 'original',
                   type = 'full', diag = T,
                   cl.cex = 1.4, tl.pos = "lt", tl.col = "black",
                   tl.cex = 1.5, tl.srt = 45,
                   # col = corrplot::COL1('YlGn', 10),
                   # is.corr = F, col.lim = c(0.8, 1),
                   addCoef.col = 'white', number.cex = 0.9,
                   mar = c(1, 1, 1, 1), title = "Pearson' correlation (Protein group level)"
)
corrplot::corrplot(cors_pep, method = 'square', order = 'original',
                   type = 'full', diag = T,
                   cl.cex = 1.4, tl.pos = "lt", tl.col = "black",
                   tl.cex = 1.5, tl.srt = 45,
                   # col = corrplot::COL1('YlGn', 10),
                   # is.corr = F, col.lim = c(0.5, 1),
                   addCoef.col = 'white', number.cex = 0.9,
                   mar = c(1, 1, 1, 1), title = "Pearson' correlation (Peptide level)"
)
graphics.off()

list(Figure2I_Protein = cors_pro %>% as.data.frame() %>% rownames_to_column(),
     Figure2I_Peptide = cors_pep %>% as.data.frame() %>% rownames_to_column()) %>% rio::export('F2I_pearson_correlation.xlsx')


## 1.11 Protein types&locations (J) --------
library(circlize)
library(ggalluvial)
my_color_palette <- function(n, ..., visible = T){
  require(RColorBrewer)
  col <- colorRampPalette(unlist(list(...)))(n)
  if(visible) {image(x = 1:n, y = 1, z = as.matrix(1:n), col = col)}
  return(col)
}
IPA_ingel <- rio::import('InGel_all_molecules.xls')
IPA_intip <- rio::import('InTip_all_molecules.xls')

colnames(IPA_ingel) <- IPA_ingel[1, ]
IPA_ingel <- IPA_ingel %>% slice(-1) %>% rename(Type = `Type(s)`)
colnames(IPA_intip) <- IPA_intip[1, ]
IPA_intip <- IPA_intip %>% slice(-1) %>% rename(Type = `Type(s)`)

### ChordDiagram ----------
# for InGel
mat_ingel <- IPA_ingel %>%
  mutate(Location = factor(Location, levels = c('Nucleus', 'Cytoplasm', 'Plasma Membrane', 'Extracellular Space', 'Other'))) %>% 
  mutate(Type = ifelse(Type == 'other', 'other types', Type)) %>% 
  reshape2::dcast(Location ~ Type) %>%
  column_to_rownames('Location') %>% 
  t()
rownames(mat_ingel) %<>% str_to_sentence()
colnames(mat_ingel) %<>% str_to_sentence()

grid.col.location <- brewer.pal(8, 'Set2')[1:5]
names(grid.col.location) <- colnames(mat_ingel)

# grid.col.type <- rep('grey', nrow(mat_ingel))
grid.col.type <- my_color_palette(nrow(mat_ingel), brewer.pal(8, 'Spectral'), visible = F)
names(grid.col.type) <- rownames(mat_ingel)
grid.col <- c(grid.col.location, grid.col.type)

# border_mat_ingel <- matrix('black', nrow = 1, ncol = ncol(mat_ingel))
# rownames(border_mat_ingel)

pdf('F2J_chordDiagram_InGel.pdf', width = 5, height = 5)
chordDiagram(mat_ingel, grid.col = grid.col,
             # annotationTrack = "grid",
             link.lwd = 2, # width
             link.lty = 2, # style
             transparency = 0.5,
             big.gap = 10,
             small.gap = 2,
             # link.border = border_mat_ingel,
)
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  # circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.7, niceFacing = TRUE)
}, bg.border = NA)
title(str_glue('IPA of InGel'), cex = 0.8)
graphics.off()
circos.clear()


# for InTip
mat_intip <- IPA_intip %>%
  mutate(Location = factor(Location, levels = c('Nucleus', 'Cytoplasm', 'Plasma Membrane', 'Extracellular Space', 'Other'))) %>% 
  mutate(Type = ifelse(Type == 'other', 'other types', Type)) %>% 
  reshape2::dcast(Location ~ Type) %>%
  column_to_rownames('Location') %>% 
  t()
rownames(mat_intip) %<>% str_to_sentence()
colnames(mat_intip) %<>% str_to_sentence()

grid.col.location <- brewer.pal(8, 'Set2')[1:5]
names(grid.col.location) <- colnames(mat_intip)

# grid.col.type <- rep('grey', nrow(mat_intip))
grid.col.type <- my_color_palette(nrow(mat_intip), brewer.pal(8, 'Spectral'), visible = F)
names(grid.col.type) <- rownames(mat_intip)
grid.col <- c(grid.col.location, grid.col.type)

# border_mat_intip <- matrix('black', nrow = 1, ncol = ncol(mat_intip))
# rownames(border_mat_intip)

pdf('F2J_chordDiagram_InTip.pdf', width = 5, height = 5)
chordDiagram(mat_intip, grid.col = grid.col,
             # annotationTrack = "grid",
             link.lwd = 2, # width
             link.lty = 2, # style
             transparency = 0.5,
             big.gap = 10,
             small.gap = 2,
             # link.border = border_mat_intip,
)
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  # circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.7, niceFacing = TRUE)
}, bg.border = NA)
title(str_glue('IPA of InTip'), cex = 0.8)
graphics.off()
circos.clear()

plot_legend1 <- data.frame(Location = colnames(mat_ingel), y = 1) %>% 
  ggplot()+
  geom_col(aes(x = Location, y = y, fill = Location), width = 1)+
  scale_fill_manual(values = grid.col.location)+
  coord_flip()+
  theme_void()
plot_legend2 <- data.frame(Type = rownames(mat_ingel), y = 1) %>% 
  ggplot()+
  geom_col(aes(x = Type, y = y, fill = Type), width = 1)+
  scale_fill_manual(values = grid.col.type)+
  coord_flip()+
  theme_void()
plot_legend <- ggpubr::ggarrange(plot_legend1, plot_legend2, ncol = 1)
ggsave('F2J_chordDiagram_legends.pdf', plot_legend, width = 8, height = 8)

pacman::p_unload(circlize)

### Sankey -----
# for in gel
dfsankey_ingel <- IPA_ingel %>%
  mutate(Type = ifelse(Type == 'other', 'other types', Type),
         Type = str_to_sentence(Type),
         Type = factor(Type, levels = sort(unique(Type), decreasing = T) %>% union('Other types', .) %>% rev())) %>% 
  mutate(Location = factor(Location, levels = c('Nucleus', 'Cytoplasm', 'Plasma Membrane', 'Extracellular Space', 'Other'))) %>%
  count(Location, Type)

plot_sankey_ingel <- ggplot(dfsankey_ingel, aes(y = n, axis1 = Location, axis2 = Type)) +
  geom_alluvium(aes(fill = Type), alpha = 1,
                width = 6/8, knot.pos = 0, reverse = T) +
  scale_fill_manual(values = grid.col.type %>% setNames(levels(dfsankey_ingel$Type))) +
  geom_stratum(alpha = 0.2, width = 6/8, reverse = T) +
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), reverse = T) +
  scale_x_continuous(breaks = 1:2, labels = c('Location', 'Type')) +
  labs(subtitle = 'IPA of InGel')+
  theme_minimal() +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )

# for in tip
dfsankey_intip <- IPA_intip %>%
  mutate(Type = ifelse(Type == 'other', 'other types', Type),
         Type = str_to_sentence(Type),
         Type = factor(Type, levels = sort(unique(Type), decreasing = T) %>% union('Other types', .) %>% rev())) %>% 
  mutate(Location = factor(Location, levels = c('Nucleus', 'Cytoplasm', 'Plasma Membrane', 'Extracellular Space', 'Other'))) %>%
  count(Location, Type)

plot_sankey_intip <- ggplot(dfsankey_intip, aes(y = n, axis1 = Location, axis2 = Type)) +
  geom_alluvium(aes(fill = Type), alpha = 1,
                width = 6/8, knot.pos = 0, reverse = T) +
  scale_fill_manual(values = grid.col.type %>% setNames(levels(dfsankey_intip$Type))) +
  geom_stratum(alpha = 0.2, width = 6/8, reverse = T) +
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), reverse = T) +
  scale_x_continuous(breaks = 1:2, labels = c('Location', 'Type')) +
  labs(subtitle = 'IPA of InTip')+
  theme_minimal() +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
pJ_sankey <- ggpubr::ggarrange(plot_sankey_ingel, plot_sankey_intip, nrow = 1)
ggsave('F2_J_Sankey.pdf', pJ_sankey, width = 12, height = 12)


## output tables 2B-2J -------
tblBtoI <- list(Figure2B = info %>% inner_join(tbl1) %>% inner_join(tbl2) %>% select(-Note),
                Figure2C = venntbl,
                Figure2D = tbl3,
                Figure2E = tbl4,
                Figure2F = info %>% inner_join(tbl5) %>% select(-Note),
                Figure2G = tbl6,
                Figure2H_Peptide = dfcv_pep,
                Figure2H_PG = dfcv_pro,
                Figure2I_Peptide = cors_pep %>% as.data.frame() %>% rownames_to_column(),
                Figure2I_PG = cors_pro %>% as.data.frame() %>% rownames_to_column(),
                Figure2J_InGel = IPA_ingel,
                Figure2J_InTip = IPA_intip)
rio::export(tblBtoI, 'F2_B-J.xlsx')

save.image('F2_B-J.RData')



# 2.L ---------------------------------------------------------------------
rm(list = ls())
gc()

## 2.1 load data ----
#peptide matrix
dfpep1 <- rio::import('//172.16.13.136/share/members/Dongzhen/proteomEx2_benchmark_result/N20230317proteomEx_v2_detection_range_DDA/combined_peptide.tsv')
dfpep1 %<>% filter(!str_detect(Protein, '^CON_'))
dfpep2 <- rio::import('//172.16.13.136/share/members/Dongzhen/proteomEx2_benchmark_result/N20230316dongz_proteomEx2_detection_range_DIA/Detection_range_report.pr_matrix.tsv')

# pep_dda <- dfpep1 %>%
#   select(`Peptide Sequence`, matches('^In\\w+?_\\d+[mu]m_\\d+ Intensity$')) %>%
#   column_to_rownames('Peptide Sequence')
# colnames(pep_dda) %<>% str_remove(' Intensity$')
pep_dda <- dfpep1 %>%
  select(`Peptide Sequence`, matches('^In\\w+?_\\d+[mu]m_\\d+ Spectral Count$')) %>%
  column_to_rownames('Peptide Sequence')
colnames(pep_dda) %<>% str_remove(' Spectral Count$')

pep_dia <- dfpep2 %>%
  select(Stripped.Sequence, matches('InTip|InGel')) %>% 
  group_by(Stripped.Sequence) %>% 
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence')
colnames(pep_dia) %<>% str_remove('^.+\\\\')

#protein matrix
dfpro1 <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx_v2_detection_range_DDA/combined_protein.tsv')
dfpro1 %<>% filter(!str_detect(Protein, '^CON_'))
dfpro2 <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230316dongz_ProteomEx2_detection_range_DIA/Detection_range_report.pg_matrix.tsv')

# pro_dda <- dfpro1 %>%
#   mutate(tmp = `Indistinguishable Proteins` %>%
#            sapply(function(x){
#              str_split(x, ', ')[[1]] %>%
#                str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
#                str_c(collapse = ', ')
#            }),
#          `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
#   select(`Protein ID`, matches('^In\\w+?_\\d+[mu]m_\\d+ Intensity$')) %>%
#   column_to_rownames('Protein ID')
# colnames(pro_dda) %<>% str_remove(' Intensity$')
pro_dda <- dfpro1 %>%
  mutate(tmp = `Indistinguishable Proteins` %>%
           sapply(function(x){
             str_split(x, ', ')[[1]] %>%
               str_remove_all('(^sp\\|)|(\\|\\w+$)') %>% 
               str_c(collapse = ', ')
           }),
         `Protein ID` = ifelse(tmp == '', `Protein ID`, str_c(`Protein ID`, ', ', tmp))) %>% 
  select(`Protein ID`, matches('^In\\w+?_\\d+[mu]m_\\d+ Total Spectral Count$')) %>%
  column_to_rownames('Protein ID')
colnames(pro_dda) %<>% str_remove(' Total Spectral Count$')

pro_dia <- dfpro2 %>%
  select(Protein.Group, matches('InTip|InGel')) %>% 
  column_to_rownames('Protein.Group')
colnames(pro_dia) %<>% str_remove('^.+\\\\')


#sample info
info <- rio::import('F2_info.xlsx', sheet = 2) %>%
  select(-Note) %>% 
  pivot_longer(cols = -(`Label in file (Diameter)`:`Label in file (Volume)`),
               values_to = 'Label in figure', values_drop_na = T) %>% 
  select(-name) %>%
  mutate(Method = str_remove(`Label in figure`, '_.+$'),
         Size = str_remove(`Label in figure`, '^.+_'))
info_dda <- data.frame(SampleName = colnames(pep_dda),
                       'Label in figure' = str_remove(colnames(pep_dda), '_\\d+$'),
                       MS = 'DDA',
                       check.names = F)
info_dia <- data.frame(SampleName = colnames(pep_dia),
                       'Label in figure' = str_extract(colnames(pep_dia), 'In\\w+?_\\d+[um]m'),
                       MS = 'DIA',
                       check.names = F)
info <- rbind(info_dda, info_dia) %>% left_join(info)
info %<>%
  mutate(Method = str_to_title(Method) %>% str_replace('In', 'In-'),
         MS = str_c(MS, '-MS'),
         Label = str_c(Method, ' ', MS))
rio::export(info, 'F2L_info.xlsx')

#regenerate dfpep and dfpro; within protein groups
dfpep_dda <- pep_dda %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpep_dia <- pep_dia %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro_dda <- pro_dda %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)
dfpro_dia <- pro_dia %>% t() %>% as.data.frame() %>% rownames_to_column('SampleName') %>% inner_join(info, .)


## 2.2 Volume and peptide/protein counts (L) --------
dfline_pep <- c(apply(pep_dda, 2, function(y) sum(y != 0)), apply(pep_dia, 2, function(y) sum(!is.na(y)))) %>% 
  as.data.frame() %>% setNames('# peptides') %>% rownames_to_column('SampleName') %>% 
  inner_join(info, .)
dfline_pro <- c(apply(pro_dda, 2, function(y) sum(y != 0)), apply(pro_dia, 2, function(y) sum(!is.na(y)))) %>% 
  as.data.frame() %>% setNames('# proteins') %>% rownames_to_column('SampleName') %>% 
  inner_join(info, .)
tblL <- dfline <- dfline_pep %>% inner_join(dfline_pro)
# d_v_fct <- max(info$`Label in file (Diameter)`) / max(info$`Label in file (Volume)`)

dfline$`Label in file (Volume)` %<>% factor(levels = sort(unique(.)))
dfline$`Label in file (Diameter)` %<>% factor(levels = sort(unique(.)))
dfline_stat <- dfline %>%
  group_by(Label, Method, `Label in file (Volume)`, `Label in file (Diameter)`) %>%
  summarise_at(vars(`# peptides`, `# proteins`), list(mean = mean, sd = sd))


# guide functions
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}

# plot
plot_pep <- ggplot(dfline_stat, aes(x = `Label in file (Volume)`, y = `# peptides_mean`, color = Label, fill = Label))+
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.9)+
  geom_errorbar(aes(ymin = `# peptides_mean` - `# peptides_sd`, ymax = `# peptides_mean` + `# peptides_sd`), width = 0.9, position = position_dodge()) +
  labs(x = 'Volume (nL)', y = '# peptides', subtitle = '')+
  scale_color_manual(values = c('In-tip DDA-MS' = '#E41A1C', 'In-gel DDA-MS' =  '#377EB8', 'In-tip DIA-MS' = '#E41A1C', 'In-gel DIA-MS' =  '#377EB8')) +
  scale_fill_manual(values = c('In-tip DDA-MS' = '#AAAAAA', 'In-gel DDA-MS' =  '#AAAAAA', 'In-tip DIA-MS' = '#EEEEEE', 'In-gel DIA-MS' =  '#EEEEEE')) +
  scale_linetype_manual(values = c('DDA-MS' = 'dashed', 'DIA-MS' = 'solid')) +
  guides(x.sec = guide_axis_label_trans(~paste(as.character(sort(unique(info$`Label in file (Diameter)`)))))) +
  scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'right',
        axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x.top = element_text(angle = 45, vjust = 0, hjust = 0))

plot_pro <- ggplot(dfline_stat, aes(x = `Label in file (Volume)`, y = `# proteins_mean`, color = Label, fill = Label))+
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.9)+
  geom_errorbar(aes(ymin = `# proteins_mean` - `# proteins_sd`, ymax = `# proteins_mean` + `# proteins_sd`), width = 0.9, position = position_dodge()) +
  labs(x = 'Volume (nL)', y = '# proteins', subtitle = '')+
  scale_color_manual(values = c('In-tip DDA-MS' = '#E41A1C', 'In-gel DDA-MS' =  '#377EB8', 'In-tip DIA-MS' = '#E41A1C', 'In-gel DIA-MS' =  '#377EB8')) +
  scale_fill_manual(values = c('In-tip DDA-MS' = '#AAAAAA', 'In-gel DDA-MS' =  '#AAAAAA', 'In-tip DIA-MS' = '#EEEEEE', 'In-gel DIA-MS' =  '#EEEEEE')) +
  scale_linetype_manual(values = c('DDA-MS' = 'dashed', 'DIA-MS' = 'solid')) +
  guides(x.sec = guide_axis_label_trans(~paste(as.character(sort(unique(info$`Label in file (Diameter)`)))))) +
  scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'right',
        axis.text.x.bottom = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x.top = element_text(angle = 45, vjust = 0, hjust = 0))

pL <- ggpubr::ggarrange(plot_pep, plot_pro, nrow = 1, common.legend = T)
ggsave('F2L_detection_range.pdf', pL, width = 11, height = 6)
list(Figure2L_raw = tblL, Figure2L_stat = dfline_stat) %>% rio::export('F2_L.xlsx')


save.image('F2L.RData')

# 3.test ------------------------------------------------------------------
## physical-chemical ------------
library(Peptides)
# Osorio, D., Rondon-Villarreal, P. & Torres, R. Peptides: A package for data mining of antimicrobial peptides. The R Journal. 7(1), 4–14 (2015).

x <- dfpep$`Peptide Sequence`[1]

pI(x, pKscale = 'EMBOSS') # isoelectic point (pI) 
charge(x, pH = 7, pKscale = 'Lehninger') # This function computes the net charge of a protein sequence based on the Henderson-Hasselbalch equation described by Moore, D. S. (1985). The net charge can be calculated at defined pH (Default = 7) using one of the 9 pKa scales availables: Bjellqvist, Dawson, EMBOSS, Lehninger, Murray, Rodwell, Sillero, Solomon or Stryer.
instaIndex(x) # The instability index is an estimation of the stability of a protein experimentally. A protein whose instability index is smaller than 40 is predicted as stable
aIndex(x) # This function calculates the Ikai (1980) aliphatic index of a protein. The aindex is defined as the relative volume occupied by aliphatic side chains (Alanine, Valine, Isoleucine, and Leucine). It may be regarded as a positive factor for the increase of thermostability of globular proteins.
crucianiProperties(x) # This function calculates the Cruciani properties of an amino-acids sequence using the scaled principal component scores that summarize a broad set of descriptors calculated based on the interaction of each amino acid residue with several chemical groups (or "probes"), such as charged ions, methyl, hydroxyl groups, and so forth.
# The computed average of Cruciani properties of all the amino acids in the corresponding peptide sequence. Each PP represent an amino-acid property as follows:
# PP1: Polarity,
# PP2: Hydrophobicity,
# PP3: H-bonding
kideraFactors(x) # The Kidera factors are a set of ten numerical descriptors of amino acids, derived using dimension reduction techniques applied to 188 physical properties of the 20 amino acids. 
# KF1: Helix/bend preference,
# KF2: Side-chain size,
# KF3: Extended structure preference,
# KF4: Hydrophobicity,
# KF5: Double-bend preference,
# KF6: Partial specific volume,
# KF7: Flat extended preference,
# KF8: Occurrence in alpha region,
# KF9: pK-C,
# KF10: Surrounding hydrophobicity
# These factors play a crucial role in various applications, such as ligand-binding site prediction and characterizing antibodies based on physicochemical properties. The computation of Kidera factors for a protein sequence involves averaging these values.
zScales(x) # Z-scales are based on physicochemical properties of the AAs including NMR data and thin-layer chromatography (TLC) data.
# The computed average of Z-scales of all the amino acids in the corresponding peptide sequence. Each Z scale represent an amino-acid property as follows:
# Z1: Lipophilicity
# Z2: Steric properties (Steric bulk/Polarizability)
# Z3: Electronic properties (Polarity / Charge)
# Z4 and Z5: They relate electronegativity, heat of formation, electrophilicity and hardness.





## maybe not clear ----------
membpos(x) # This function calculates the theoretical class of a protein sequence based on the relationship between the hydrophobic moment and hydrophobicity scale proposed by Eisenberg (1984).
mswhimScores(x) # MS-WHIM indexes are a collection of 36 statistical indexes aimed at extracting and condensing steric and electrostatic 3D-properties of a molecule.
stScales(x) # ST-scales were proposed by Yang et al, taking 827 properties into account which are mainly constitutional, topological, geometrical, hydrophobic, elec- tronic, and steric properties of a total set of 167 AAs.
tScales(x) # T-scales are based on 67 common topological descriptors of 135 amino acids. These topological descriptors are based on the connectivity table of amino acids alone, and to not explicitly consider 3D properties of each structure.
vhseScales(x) # VHSE-scales (principal components score Vectors of Hydrophobic, Steric, and Electronic properties), is derived from principal components analysis (PCA) on independent families of 18 hydrophobic properties, 17 steric properties, and 15 electronic properties, respectively, which are included in total 50 physicochemical variables of 20 coded amino acids.
protFP(x) # The ProtFP descriptor set was constructed from a large initial selection of indices obtained from the AAindex database for all 20 naturally occurring amino acids.


aaDescriptors(x)
# The function return 66 amino acid descriptors for the 20 natural amino acids. Available descriptors are:
#   
# crucianiProperties: Cruciani, G., Baroni, M., Carosati, E., Clementi, M., Valigi, R., and Clementi, S. (2004) Peptide studies by means of principal properties of amino acids derived from MIF descriptors. J. Chemom. 18, 146-155.,
# 
# kideraFactors: Kidera, A., Konishi, Y., Oka, M., Ooi, T., & Scheraga, H. A. (1985). Statistical analysis of the physical properties of the 20 naturally occurring amino acids. Journal of Protein Chemistry, 4(1), 23-55.,
# 
# zScales: Sandberg M, Eriksson L, Jonsson J, Sjostrom M, Wold S: New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids. J Med Chem 1998, 41:2481-2491.,
# 
# FASGAI: Liang, G., & Li, Z. (2007). Factor analysis scale of generalized amino acid information as the source of a new set of descriptors for elucidating the structure and activity relationships of cationic antimicrobial peptides. Molecular Informatics, 26(6), 754-763.,
# 
# tScales: Tian F, Zhou P, Li Z: T-scale as a novel vector of topological descriptors for amino acids and its application in QSARs of peptides. J Mol Struct. 2007, 830: 106-115. 10.1016/j.molstruc.2006.07.004.,
# 
# VHSE: VHSE-scales (principal components score Vectors of Hydrophobic, Steric, and Electronic properties), is derived from principal components analysis (PCA) on independent families of 18 hydrophobic properties, 17 steric properties, and 15 electronic properties, respectively, which are included in total 50 physicochemical variables of 20 coded amino acids.,
# 
# protFP: van Westen, G. J., Swier, R. F., Wegner, J. K., IJzerman, A. P., van Vlijmen, H. W., & Bender, A. (2013). Benchmarking of protein descriptor sets in proteochemometric modeling (part 1): comparative study of 13 amino acid descriptor sets. Journal of cheminformatics, 5(1), 41.,
# 
# stScales: Yang, L., Shu, M., Ma, K., Mei, H., Jiang, Y., & Li, Z. (2010). ST-scale as a novel amino acid descriptor and its application in QSAM of peptides and analogues. Amino acids, 38(3), 805-816.,
# 
# BLOSUM: Georgiev, A. G. (2009). Interpretable numerical descriptors of amino acid space. Journal of Computational Biology, 16(5), 703-723.,
# 
# MSWHIM: Zaliani, A., & Gancia, E. (1999). MS-WHIM scores for amino acids: a new 3D-description for peptide QSAR and QSPR studies. Journal of chemical information and computer sciences, 39(3), 525-533.


aaComp(x) # amino-acid composition
# The output is a matrix with the number and percentage of amino acids of a particular class:
# Tiny (A + C + G + S + T)
# Small (A + B + C + D + G + N + P + S + T + V)
# Aliphatic (A + I + L + V)
# Aromatic (F + H + W + Y)
# Non-polar (A + C + F + G + I + L + M + P + V + W + Y)
# Polar (D + E + H + K + N + Q + R + S + T + Z)
# Charged (B + D + E + H + K + R + Z)
# Basic (H + K + R)
# Acidic (B + D + E + Z)

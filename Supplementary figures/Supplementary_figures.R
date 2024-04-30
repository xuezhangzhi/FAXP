# 0.Set environment ----------
#rm(list = ls())
#gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)


# 1.SFigure1 --------------------------------------------------------------
my_proj <- c('ProteomEx (5.6 nL)', 'Current (1.55 nL)')

df1 <- rio::import('data/Supp Fig1.xlsx')
tblS1 <- df1[c(1, 2, 3, 6)] %>%
  setNames(c('Project', 'Organism', '# peptides', '# proteins'))
tblS1[tblS1$Project == 'v1', 'Project'] <- my_proj[1]
tblS1[tblS1$Project == 'v2', 'Project'] <- my_proj[2]
tblS1$Project %<>% factor(levels = my_proj)
dfbar <- tblS1 %>% pivot_longer(cols = -(Project:Organism), names_to = 'Type', values_to = 'Number')

plotS1 <- ggplot(dfbar, aes(x = Project, y = Number, color = Project))+
  facet_wrap('Type', scales = 'free_y'
             ) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = list(my_proj),
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF1.pdf', plotS1, width = 6, height = 6)


# 2.SFigure3 --------------------------------------------------------------
#read data
info2 <- rio::import('data/Supp Fig3.xlsx')
pep_ls <- list.files('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20211215dongz_ProteomEx2.0_FASPv7_result_variable_DDA', '^peptide\\.tsv$', full.names = T, recursive = T)
pro_ls <- list.files('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20211215dongz_ProteomEx2.0_FASPv7_result_variable_DDA', '^protein\\.tsv$', full.names = T, recursive = T)
pep_ls %<>% str_subset(str_c(info2$`Name in the folder`, collapse = '|'))
pro_ls %<>% str_subset(str_c(info2$`Name in the folder`, collapse = '|'))
df2pep_ls <- lapply(pep_ls, rio::import)
df2pro_ls <- lapply(pro_ls, rio::import)

fileinfo2 <- data.frame(SampleName = c(pep_ls, pro_ls), Type = c(rep('pep', length(pep_ls)), rep('pro', length(pro_ls))))
fileinfo2$ID <- str_extract(fileinfo2$SampleName, str_c(info2$`Name in the folder`, '_\\d+', collapse = '|'))
fileinfo2$`Name in the folder` <- str_extract(fileinfo2$ID, str_c(info2$`Name in the folder`, collapse = '|'))
fileinfo2 %<>% left_join(info2)
names(df2pep_ls) <- fileinfo2$ID[fileinfo2$Type == 'pep']
names(df2pro_ls) <- fileinfo2$ID[fileinfo2$Type == 'pro']

#peptide matrix
df2pep <- plyr::ldply(df2pep_ls, .id = 'ID')
df2pep %<>% filter(!str_detect(Protein, '^CON_'))
df2pep %<>% filter(!str_detect(`Mapped Proteins`, '^CON_'))
mat2pep_count <- df2pep %>% pivot_wider(id_cols = ID, names_from = Peptide, values_from = `Spectral Count`) %>% column_to_rownames('ID') %>% t()

#protein matrix
df2pro <- plyr::ldply(df2pro_ls, .id = 'ID')
df2pro %<>% filter(!str_detect(Protein, '^CON_'))
df2pro %<>% filter(!str_detect(`Indistinguishable Proteins`, '^CON_'))
mat2pro_count <- df2pro %>% pivot_wider(id_cols = ID, names_from = `Protein ID`, values_from = `Total Spectral Count`) %>% column_to_rownames('ID') %>% t()


## 2.1 numbers (B) -----
# # peptides, # proteins
# % cysteine contained peptides, % carbamidomethyl mod
npep <- apply(mat2pep_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# peptides') %>% rownames_to_column('ID')
npro <- apply(mat2pro_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# proteins') %>% rownames_to_column('ID')
nC <- df2pep %>%
  group_by(ID) %>%
  summarise(`% cysteine contained peptides` = 100 * sum(str_detect(Peptide, 'C')) / length(Peptide))
nCmod <- df2pep %>% filter(str_detect(Peptide, 'C')) %>%
  group_by(ID) %>%
  summarise(`% carbamidomethyl modification` = 100 * sum(str_detect(`Assigned Modifications`, 'C')) / length(`Assigned Modifications`))


tblS3B <- fileinfo2 %>% select(-(SampleName:Type)) %>% distinct() %>% 
  left_join(npep) %>% left_join(npro) %>% left_join(nC) %>% left_join(nCmod)
tblS3B$`Name used in the image` %<>% factor(levels = info2$`Name used in the image`)
tblS3B %<>% arrange(`Name used in the image`)
dfbar <- tblS3B %>%
  pivot_longer(cols = -(ID:`Name used in the image`), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2), c(1,3), c(2,3)), function(indice){
  levels(tblS3B$`Name used in the image`)[indice]
})
plotS3B <- ggplot(dfbar, aes(x = `Name used in the image`, y = Number, color = `Name used in the image`))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF3B_pep_pro_C_Cmod.pdf', plotS3B, width = 10, height = 7)

## 2.2 % missed cleavage (C) -----
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

df_missCleav <- calc_missCleav(df2pep$Peptide)

tblS3C <- df2pep %>% left_join(df_missCleav %>% rename(Peptide = pepseq), relationship = 'many-to-many') %>%
  count(ID, MissedCleavage) %>%
  with_groups(ID, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  inner_join(tblS3B %>% select(ID:`Name used in the image`), .)

plotS3C <- ggplot(tblS3C, aes(x = `Name used in the image`, y = `ratio (%)`, color = `Name used in the image`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free_y',) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF')) +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF3C_missed_cleavage.pdf', plotS3C, width = 10, height = 7)


## 2.3 peptide length statistics (D) -----
tblS3D <- df2pep %>%
  inner_join(fileinfo2) %>% 
  distinct(`Name used in the image`, Peptide, `Peptide Length`) %>% 
  count(`Name used in the image`, `Peptide Length`)

plotS3D <- ggplot(tblS3D, aes(x = `Peptide Length`, y = n, color = `Name used in the image`, fill = `Name used in the image`)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = '', y = '', subtitle = 'Length of peptides')+
  # scale_y_continuous(labels = scales::scientific) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')
ggsave('SF3D_peptide_length.pdf', plotS3D, width = 5, height = 5)


# 3.SFigure4 --------------------------------------------------------------
info3 <- rio::import('data/Supp Fig4.xlsx')
pep_ls <- list.files('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20211215dongz_ProteomEx2.0_FASPv7_result_fixed_DDA', '^peptide\\.tsv$', full.names = T, recursive = T)
pro_ls <- list.files('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20211215dongz_ProteomEx2.0_FASPv7_result_fixed_DDA', '^protein\\.tsv$', full.names = T, recursive = T)
pep_ls %<>% str_subset(str_c(info3$`Name in the folder`, collapse = '|'))
pro_ls %<>% str_subset(str_c(info3$`Name in the folder`, collapse = '|'))
df3pep_ls <- lapply(pep_ls, rio::import)
df3pro_ls <- lapply(pro_ls, rio::import)

fileinfo3 <- data.frame(SampleName = c(pep_ls, pro_ls), Type = c(rep('pep', length(pep_ls)), rep('pro', length(pro_ls))))
fileinfo3$ID <- str_extract(fileinfo3$SampleName, str_c(info3$`Name in the folder`, '_\\d+', collapse = '|'))
fileinfo3$`Name in the folder` <- str_extract(fileinfo3$ID, str_c(info3$`Name in the folder`, collapse = '|'))
fileinfo3 %<>% left_join(info3)
names(df3pep_ls) <- fileinfo3$ID[fileinfo3$Type == 'pep']
names(df3pro_ls) <- fileinfo3$ID[fileinfo3$Type == 'pro']

#peptide matrix
df3pep <- plyr::ldply(df3pep_ls, .id = 'ID')
df3pep %<>% filter(!str_detect(Protein, '^CON_'))
df3pep %<>% filter(!str_detect(`Mapped Proteins`, '^CON_'))
mat3pep_count <- df3pep %>% pivot_wider(id_cols = ID, names_from = Peptide, values_from = `Spectral Count`) %>% column_to_rownames('ID') %>% t()

#protein matrix
df3pro <- plyr::ldply(df3pro_ls, .id = 'ID')
df3pro %<>% filter(!str_detect(Protein, '^CON_'))
df3pro %<>% filter(!str_detect(`Indistinguishable Proteins`, '^CON_'))
mat3pro_count <- df3pro %>% pivot_wider(id_cols = ID, names_from = `Protein ID`, values_from = `Total Spectral Count`) %>% column_to_rownames('ID') %>% t()


## 3.1 numbers (B) -----
# # peptides, # proteins
# % cysteine contained peptides, % carbamidomethyl mod
npep <- apply(mat3pep_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# peptides') %>% rownames_to_column('ID')
npro <- apply(mat3pro_count, 2, function(y) sum(!is.na(y))) %>% 
  as.data.frame() %>% setNames('# proteins') %>% rownames_to_column('ID')
nC <- df3pep %>%
  group_by(ID) %>%
  summarise(`% cysteine contained peptides` = 100 * sum(str_detect(Peptide, 'C')) / length(Peptide))
nCmod <- df3pep %>% filter(str_detect(Peptide, 'C')) %>%
  group_by(ID) %>%
  summarise(`% carbamidomethyl modification` = 100 * sum(str_detect(`Assigned Modifications`, 'C')) / length(`Assigned Modifications`))


tblS4B <- fileinfo3 %>% select(-(SampleName:Type)) %>% distinct() %>% 
  left_join(npep) %>% left_join(npro) %>% left_join(nC) %>% left_join(nCmod)
tblS4B$`Name used in the image` %<>% factor(levels = info3$`Name used in the image`)
tblS4B %<>% arrange(`Name used in the image`)
dfbar <- tblS4B %>%
  pivot_longer(cols = -(ID:`Name used in the image`), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2), c(1,3), c(2,3)), function(indice){
  levels(tblS4B$`Name used in the image`)[indice]
})
plotS4B <- ggplot(dfbar, aes(x = `Name used in the image`, y = Number, color = `Name used in the image`))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF4B_pep_pro_C_Cmod.pdf', plotS4B, width = 10, height = 5)

## 3.2 % missed cleavage (C) -----
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

df_missCleav <- calc_missCleav(df3pep$Peptide)

tblS4C <- df3pep %>% left_join(df_missCleav %>% rename(Peptide = pepseq), relationship = 'many-to-many') %>%
  count(ID, MissedCleavage) %>%
  with_groups(ID, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  inner_join(tblS4B %>% select(ID:`Name used in the image`), .)

plotS4C <- ggplot(tblS4C, aes(x = `Name used in the image`, y = `ratio (%)`, color = `Name used in the image`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free_y',) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF')) +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF4C_missed_cleavage.pdf', plotS4C, width = 10, height = 7)


## 3.3 peptide length statistics (D) -----
tblS4D <- df3pep %>%
  inner_join(fileinfo3) %>% 
  distinct(`Name used in the image`, Peptide, `Peptide Length`) %>% 
  count(`Name used in the image`, `Peptide Length`)

plotS4D <- ggplot(tblS4D, aes(x = `Peptide Length`, y = n, color = `Name used in the image`, fill = `Name used in the image`)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = '', y = '', subtitle = 'Length of peptides')+
  # scale_y_continuous(labels = scales::scientific) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')
ggsave('SF4D_peptide_length.pdf', plotS4D, width = 5, height = 5)

# 4.SFigure6 --------------------------------------------------------------
## 4.1 read data ------
info4 <- rio::import('data/20230219_Batch_design_CRC_demo.xlsx', sheet = 'Sheet2') %>% select(1:7)
df4pg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_demo_result/N20230320ProteomEx_v2_demo_whole_library_Yueliang/Whole_lib_report.pg_matrix.tsv')
df4pr <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_demo_result/N20230320ProteomEx_v2_demo_whole_library_Yueliang/Whole_lib_report.pr_matrix.tsv')

colnames(df4pg) %<>% str_remove('^.+\\\\')
colnames(df4pr) %<>% str_remove('^.+\\\\')

info4 <- data.frame(
  SampleName = colnames(df4pg)[-(1:5)],
  Batch = str_extract(colnames(df4pg)[-(1:5)], 'b\\d+_(\\d+|pool)')
) %>% full_join(info4)

mat4pg <- df4pg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat4pep <- df4pr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

## 4.2 pool QC ------
### 4.2.1 protein, peptide count (B)-----
df4_count <- data.frame(`# proteins` = apply(mat4pg, 2, function(y) sum(!is.na(y))),
                     `# peptides` = apply(mat4pep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info4, .)

tblS6B <- dfbar <- df4_count %>%
  filter(str_detect(Batch, 'pool')) %>% 
  select(SampleName, Batch, `# proteins`, `# peptides`) %>% 
  pivot_longer(cols = -(SampleName:Batch), names_to = 'Type', values_to = 'Number') %>% 
  mutate(Label = 'Pool')

plotS6B <- ggplot(dfbar, aes(x = Label, y = Number))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none')
ggsave('SF6B_pep_pro_number.pdf', plotS6B, width = 4, height = 5)

### 4.2.2 CV (C) ---------
tblS6C_protein <- df4pg_poolcv <- mat4pg %>% t() %>% as.data.frame() %>%
  rownames_to_column('SampleName') %>%
  inner_join(info4, .) %>%
  filter(str_detect(Batch, 'pool')) %>%
  select(-(ID:Random), -SampleName) %>% 
  pivot_longer(cols = -Batch, names_to = 'Protein', values_to = 'Intensity', values_drop_na = T) %>% 
  group_by(Protein) %>% 
  mutate(Intensity = 2^Intensity) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean, Label = 'Pool')

plotS6C_pro <- ggplot(df4pg_poolcv, aes(x = Label, y = CV))+
  geom_violin(fill = '#AAAAAA', color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Protein group')+
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

tblS6C_peptide <- df4pep_poolcv <- mat4pep %>% t() %>% as.data.frame() %>%
  rownames_to_column('SampleName') %>%
  inner_join(info4, .) %>%
  filter(str_detect(Batch, 'pool')) %>%
  select(-(ID:Random), -SampleName) %>% 
  pivot_longer(cols = -Batch, names_to = 'Peptide', values_to = 'Intensity', values_drop_na = T) %>% 
  group_by(Peptide) %>% 
  mutate(Intensity = 2^Intensity) %>% 
  summarise(mean = mean(Intensity),
            sd = sd(Intensity)) %>% 
  mutate(CV = sd / mean, Label = 'Pool')

plotS6C_pep <- ggplot(df4pep_poolcv, aes(x = Label, y = CV))+
  geom_violin(fill = '#AAAAAA', color = NA, width = 1)+
  # stat_boxplot(geom = 'errorbar', color = '#000000', width = 0.05)+
  geom_boxplot(fill = '#FFFFFF', color = '#000000', width = 0.1, outlier.size = 0.5)+
  labs(x = '', y = 'CV', subtitle = 'Peptide')+
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none')

plotS6C <- ggpubr::ggarrange(plotS6C_pep, plotS6C_pro, nrow = 1, ncol = 2)
ggsave('SF6C_pep_pro_CV.pdf', plotS6C, width = 5, height = 5)


### 4.2.3 Pearson's correlation (D) ---------
info4_pool <- info4 %>% filter(str_detect(Batch, 'pool'))

cors_pro <- mat4pg %>%
  select(all_of(info4_pool$SampleName)) %>%
  set_colnames(info4_pool$Batch) %>% 
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>%
  round(2)
min(cors_pro) # 0.97

cors_pep <- mat4pep %>%
  select(all_of(info4_pool$SampleName)) %>%
  set_colnames(info4_pool$Batch) %>% 
  cor(use = 'pairwise.complete.obs', method = 'pearson') %>%
  round(2)
min(cors_pep) # 0.93

pdf('SF6D_pearson_correlation.pdf', width = 5, height = 5)
corrplot::corrplot(cors_pro, method = 'square', order = 'original',
                   type = 'full', diag = T,
                   cl.cex = 1.4, tl.pos = "lt", tl.col = "black",
                   tl.cex = 1.5, tl.srt = 45,
                   # col = corrplot::COL1('YlGn', 10),
                   # is.corr = F, col.lim = c(0.8, 1),
                   addCoef.col = 'white', number.cex = 0.9,
                   mar = c(1, 1, 1, 1), title = "Pearson's correlation (Protein group level)"
)
corrplot::corrplot(cors_pep, method = 'square', order = 'original',
                   type = 'full', diag = T,
                   cl.cex = 1.4, tl.pos = "lt", tl.col = "black",
                   tl.cex = 1.5, tl.srt = 45,
                   # col = corrplot::COL1('YlGn', 10),
                   # is.corr = F, col.lim = c(0.5, 1),
                   addCoef.col = 'white', number.cex = 0.9,
                   mar = c(1, 1, 1, 1), title = "Pearson's correlation (Peptide level)"
)
graphics.off()


## 4.3 missed cleavage (DDA data) (E) ----
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

ddapeps4 <- list.files('//172.16.13.136/share/members/Dongzhen/ProteomEx2_demo_result/N20230322ProteomEx_v2_CRC_DDA', pattern = '^peptide\\.tsv$', recursive = T, full.names = T)
names(ddapeps4) <- str_extract(ddapeps4, 'exp_\\d+')

df4pepdda <- plyr::ldply(ddapeps4, rio::import, .id = 'ID')
# df4pepdda %>% filter(str_detect(Protein, 'CON_'))
# df4pepdda %>% filter(str_detect(`Mapped Proteins`, 'CON_'))
pepList <- df4pepdda %>% count(ID, Peptide) %>% 
  plyr::dlply('ID', function(dfsub){
    dfsub$Peptide
  })
df_missCleav <- plyr::ldply(pepList, calc_missCleav, .id = 'SampleName')
tblS6E <- df_missCleav %>%
  count(ID, MissedCleavage) %>%
  with_groups(ID, mutate, `ratio (%)` = 100 * n / sum(n)) %>% 
  mutate(Label = 'ProteomEx_v2_CRC_DDA')

plotS6E <- ggplot(tblS6E, aes(x = Label, y = `ratio (%)`)) +
  facet_wrap(vars(str_c('Missed cleavage=', MissedCleavage)), scales = 'free_y',) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '% missed cleavage peptides')+
  theme_classic() +
  theme(text = element_text(size = 15), legend.position = 'none',
        strip.background = element_rect(fill = '#FFFFFF'))
ggsave('SF6E_missed_cleavage.pdf', plotS6E, width = 8, height = 4)


## 4.4 quantified proteins (DIA data) (F) ----
df4pg_sample <- mat4pg %>% t() %>% as.data.frame() %>%
  rownames_to_column('SampleName') %>%
  inner_join(info4, .) %>%
  filter(!str_detect(Batch, 'pool'))

proList_patient <- plyr::dlply(df4pg_sample, 'Patient', function(dfsub){
  mat_tmp <- dfsub %>%
    column_to_rownames('Batch') %>%
    select(-(SampleName:Random)) %>% 
    t() %>% as.data.frame()
  mat_tmp <- mat_tmp[apply(mat_tmp, 1, function(x) any(!is.na(x))), ]
  return(rownames(mat_tmp))
})

proList_slide <- df4pg_sample %>% #of patient1
  filter(Patient == 'P1') %>% 
  plyr::dlply('Slide', function(dfsub){
  mat_tmp <- dfsub %>%
    column_to_rownames('Batch') %>%
    select(-(SampleName:Random)) %>% 
    t() %>% as.data.frame()
  mat_tmp <- mat_tmp[apply(mat_tmp, 1, function(x) any(!is.na(x))), ]
  return(rownames(mat_tmp))
})
names(proList_slide) %<>% c('1' = 'P1S1', '2' = 'P1S2', '3' = 'P1S3')[.]

proList_region <- df4pg_sample %>% #of patient1 slide1
  filter(Patient == 'P1', Slide == 1) %>% 
  mutate(Region = factor(Region, levels = c('N', 'L', 'H', 'PC', 'CC'))) %>% 
  plyr::dlply('Region', function(dfsub){
    mat_tmp <- dfsub %>%
      column_to_rownames('Batch') %>%
      select(-(SampleName:Random)) %>% 
      t() %>% as.data.frame()
    mat_tmp <- mat_tmp[apply(mat_tmp, 1, function(x) any(!is.na(x))), ]
    return(rownames(mat_tmp))
  })
names(proList_region) %<>% str_c('P1S1', .)

pdf('SF6F_upset1.pdf', width = 5, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(proList_patient),
    sets = rev(names(proList_patient)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on', order.by = 'degree', # order.by = 'freq',
    mainbar.y.label = '# proteins',
    point.size = 3, line.size = 1, set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5), mb.ratio = c(0.6, 0.4)
  ))
print(
  UpSetR::upset(
    UpSetR::fromList(proList_slide),
    sets = rev(names(proList_slide)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on', order.by = 'degree', # order.by = 'freq',
    mainbar.y.label = '# proteins',
    point.size = 3, line.size = 1, set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5), mb.ratio = c(0.6, 0.4)
  ))
graphics.off()

pdf('SF6F_upset2.pdf', width = 10, height = 6)
print(
  UpSetR::upset(
    UpSetR::fromList(proList_region),
    sets = rev(names(proList_region)), keep.order = T, # from UpSet bottom to top
    empty.intersections = 'on', order.by = 'degree', # order.by = 'freq',
    mainbar.y.label = '# proteins',
    point.size = 3, line.size = 1, set_size.show = T, number.angles = 0,
    text.scale = c(2, 2, 2, 2, 2, 1.5), mb.ratio = c(0.6, 0.4)
  ))
graphics.off()


list2dfcoerce <- function(mylist){
  # should be a flatten list
  n <- length(mylist)
  m <- max(sapply(mylist, length))
  res <- matrix(nrow = m, ncol = n)
  for(j in 1:n){
    for(i in 1:m){
      res[i, j] <- mylist[[j]][i]
    }
  }
  res <- as.data.frame(res)
  colnames(res) <- names(mylist)
  # variable_class <- apply(res, 2, class)
  # funcs <- sapply(paste0('as.', variable_class), match.fun)
  # res[] <- Map(function(dd, f) f(as.character(dd)), res, funcs)
  return(res)
}
tbls6F <- proList_patient %>%
  append(proList_slide) %>%
  append(proList_region) %>% 
  list2dfcoerce()


# 5.SFigure7 --------------------------------------------------------------
## 5.1 Fig.S7A gel-making ------
df5Apg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317dongz_ProteomEx2_homogenization_condition_DIA/Homogenization_report.pg_matrix.tsv')
df5Apr <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317dongz_ProteomEx2_homogenization_condition_DIA/Homogenization_report.pr_matrix.tsv')
colnames(df5Apg) %<>% str_remove('^.+\\\\')
colnames(df5Apr) %<>% str_remove('^.+\\\\')

info5A <- data.frame(SampleName = colnames(df5Apg)[-(1:5)])
info5A$Label <- str_extract(info5A$SampleName, 'Gel\\d+')

mat5Apg <- df5Apg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat5Apep <- df5Apr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

tblS7A <- data.frame(`# proteins` = apply(mat5Apg, 2, function(y) sum(!is.na(y))),
                     `# peptides` = apply(mat5Apep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info5A, .)

dfbar <- tblS7A %>%
  pivot_longer(cols = -(SampleName:Label), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)), function(indice){
  unique(tblS7A$Label)[indice]
})
plotS7A <- ggplot(dfbar, aes(x = Label, y = Number, color = Label))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF7A_pep_pro.pdf', plotS7A, width = 7, height = 5)


## 5.2 Fig.S7B LC gradient ------
df5Bpg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230316dongz_ProteomEx2_LC_gradients_DIA/LC_gradients_report.pg_matrix.tsv')
df5Bpr <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230316dongz_ProteomEx2_LC_gradients_DIA/LC_gradients_report.pr_matrix.tsv')
colnames(df5Bpg) %<>% str_remove('^.+\\\\')
colnames(df5Bpr) %<>% str_remove('^.+\\\\')

info5B <- rio::import('data/Supp Fig7.xlsx')
info5B <- data.frame(
  SampleName = colnames(df5Bpg)[-(1:5)],
  `Name in the folder` = str_extract(colnames(df5Bpg)[-(1:5)], str_c(info5B$`Name in the folder`, collapse = '|')),
  check.names = F
) %>% inner_join(info5B)

mat5Bpg <- df5Bpg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat5Bpep <- df5Bpr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

tblS7B <- data.frame(`# proteins` = apply(mat5Bpg, 2, function(y) sum(!is.na(y))),
                     `# peptides` = apply(mat5Bpep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info5B, .)

# add FWHM.Scans
df5Bstat <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230316dongz_ProteomEx2_LC_gradients_DIA/LC_gradients_report.stats.tsv')
df5Bstat$File.Name %<>% str_remove('^.+\\\\')
tblS7B <- df5Bstat %>% select(File.Name, FWHM.Scans) %>%
  rename(SampleName = File.Name) %>% 
  inner_join(tblS7B, .)

dfbar <- tblS7B %>%
  pivot_longer(cols = -(SampleName:`Name used in the image`), names_to = 'Type', values_to = 'Number')

my_comparisons <- lapply(list(c(1,2)), function(indice){
  unique(tblS7B$`Name used in the image`)[indice]
})
plotS7B <- ggplot(dfbar, aes(x = `Name used in the image`, y = Number, color = `Name used in the image`))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none') +
  ggpubr::stat_compare_means(
    method = 't.test',
    map_signif_level = F,# hide.ns = T, label = 'p.signif',
    comparisons = my_comparisons,
    size = 5, hjust = 0.5, vjust = 0)
ggsave('SF7B_pep_pro_FWHM.pdf', plotS7B, width = 7.5, height = 5)




# 6.SFigure8 --------------------------------------------------------------
df6pg <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx2_minute_sample_DIA/minute_sample_report.pg_matrix.tsv')
df6pr <- rio::import('//172.16.13.136/share/members/Dongzhen/ProteomEx2_benchmark_result/N20230317ProteomEx2_minute_sample_DIA/minute_sample_report.pr_matrix.tsv')
colnames(df6pg) %<>% str_remove('^.+\\\\')
colnames(df6pr) %<>% str_remove('^.+\\\\')

info6 <- data.frame(SampleName = colnames(df6pg)[-(1:5)])
info6$Label <- str_extract(info6$SampleName, 'InTip_10cells')

mat6pg <- df6pg %>%
  column_to_rownames('Protein.Ids') %>% 
  select(-(Protein.Group:First.Protein.Description)) %>%
  log2()

mat6pep <- df6pr %>%
  group_by(Stripped.Sequence) %>% 
  select(-(Protein.Group:Precursor.Id)) %>%
  summarise_all(mean, na.rm = T) %>% 
  column_to_rownames('Stripped.Sequence') %>%
  log2()

tblS8B <- data.frame(`# proteins` = apply(mat6pg, 2, function(y) sum(!is.na(y))),
                     `# peptides` = apply(mat6pep, 2, function(y) sum(!is.na(y))), check.names = F) %>% 
  rownames_to_column('SampleName') %>% 
  inner_join(info6, .)

dfbar <- tblS8B %>%
  pivot_longer(cols = -(SampleName:Label), names_to = 'Type', values_to = 'Number')

plotS8B <- ggplot(dfbar, aes(x = Label, y = Number, color = Label))+
  facet_wrap('Type', scales = 'free_y', nrow = 1) +
  stat_boxplot(geom = 'errorbar', width = 0.5)+
  geom_boxplot(fill = '#FFFFFF', width = 0.6, outlier.shape = NA)+
  geom_jitter(alpha = 1, size = 3, width = 0.35)+
  labs(x = '', y = '', subtitle = '')+
  scale_color_brewer(palette = 'Set1') +
  # scale_y_continuous(labels = scales::scientific) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'none')
ggsave('SF8B_pep_pro.pdf', plotS8B, width = 4, height = 5)





# Output ------------------------------------------------------------------
list(
  SFigure1 = tblS1,
  SFigure3B = tblS3B,
  SFigure3C = tblS3C,
  SFigure3D = tblS3D,
  SFigure4B = tblS4B,
  SFigure4C = tblS4C,
  SFigure4D = tblS4D,
  SFigure6B = tblS6B,
  SFigure6C_peptide = tblS6C_peptide,
  SFigure6C_protein = tblS6C_protein,
  SFigure6D_peptide = cors_pep %>% as.data.frame() %>% rownames_to_column(),
  SFigure6D_protein = cors_pro %>% as.data.frame() %>% rownames_to_column(),
  SFigure6E = tblS6E,
  SFigure6F = tbls6F,
  SFigure7A = tblS7A,
  SFigure7B = tblS7B,
  SFigure8B = tblS8B
) %>%
  rio::export('Supplementary_figures.xlsx')

save.image('Supplementary_figures.RData')




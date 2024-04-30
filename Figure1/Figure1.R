# 0.Set environment ----------
#rm(list = ls())
#gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(showtext)



# 1.read data -------------------------------------------------------------
df <- rio::import('20231129Figure6B_timeline.xlsx')
df_proteomex <- df[1:10, 1:5] %>%
  set_colnames(.[1, ]) %>% slice(-1) %>% 
  # mutate(ProteomEx = factor(ProteomEx, levels = ProteomEx)) %>% 
  mutate_at(vars(`T1 (hours)`:`Time lasting (hours)`), as.double)
df_faxp <- df[1:12, 7:11] %>%
  set_colnames(.[1, ]) %>% slice(-1) %>% 
  # mutate(FAXP = factor(FAXP, levels = FAXP)) %>% 
  mutate_at(vars(`T1 (hours)`:`Time lasting (hours)`), as.double)
progress_all <- union(df_proteomex$ProteomEx, df_faxp$FAXP) %>% 
  factor(levels = c('Sectioning & de-waxing', 'Anchoring', 'Monomer incubation', 'Gelation', 'Homogenization', 'Reduction & alkylation (whole gel)', 'Staining & expansion', 'Imaging', 'Micro-dissection', 'Reduction and alkylation (punch), in-gel digestion & peptide extraction', 'Filter-aided in-gel digestion & peptide elution', 'LC-MS/MS')) %>% sort()
proteomex_only <- setdiff(progress_all, df_faxp$FAXP)
faxp_only <- setdiff(progress_all, df_proteomex$ProteomEx)


dfbar <- list(ProteomEx = df_proteomex %>% rename(Process = ProteomEx),
              FAXP = df_faxp %>% rename(Process = FAXP)) %>% 
  plyr::ldply(.id = 'Project') %>% 
  mutate(Process = factor(Process, levels = rev(levels(progress_all))))
dfbar[dfbar$Process %in% proteomex_only, 'Label'] <- 'ProteomEx only'
dfbar[dfbar$Process %in% faxp_only, 'Label'] <- 'FAXP only'
# dfbar[is.na(dfbar$Label), 'Label'] <- 'NA'
dfbar %<>% arrange(Project, desc(Process))

# set.seed(10)
# fill_colors <- sample(brewer.pal(12, 'Paired'))
plot6 <- ggplot(dfbar, aes(x = Project, y = `Time lasting (hours)`, color = Label, group = Process, fill = Process)) +
  geom_bar(position = 'stack', stat = 'identity', width = 0.9, linewidth = 0.5) +
  labs(x = '', y = 'Timeline (hours)', subtitle = '')+
  coord_flip() +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_color_manual(values = c('ProteomEx only' = 'black', 'FAXP only' = 'black'), na.value = '#DDDDDD') +
  scale_fill_brewer(palette = 'Paired') +
  # scale_fill_manual(values = fill_colors) +
  theme_classic() +
  theme(text = element_text(size = 10))
ggsave('F6B_timeline.pdf', plot6, width = 12, height = 3)
list(Figure6B = dfbar) %>% rio::export('F6B_timeline.xlsx')


# 2.new style -2023.12.26 -----------
num2chr <- c('①', '②', '③', '④', '⑤', '⑥', '⑦', '⑧', '⑨', '⑩', '⑪', '⑫') %>%
  setNames(1:12)
dfplot <- list(ProteomEx = df_proteomex %>% rename(Process = ProteomEx),
               FAXP = df_faxp %>% rename(Process = FAXP)) %>% 
  plyr::ldply(.id = 'Project') %>% 
  mutate(Process = factor(Process, levels = levels(progress_all))) %>% 
  add_column(x = as.numeric(.$Process) %>% num2chr[.], .before = 1)
dfplot[dfplot$Process %in% proteomex_only, 'Label'] <- 'ProteomEx only'
dfplot[dfplot$Process %in% faxp_only, 'Label'] <- 'FAXP only'
dfplot %<>% arrange(Project, x)

dfplot %<>% group_by(x) %>%
  reframe(Project = Project, y2 = `Time lasting (hours)` / sum(`Time lasting (hours)`)) %>%
  ungroup() %>% 
  left_join(dfplot, .)

p1 <- ggplot(dfplot, aes(x = x, y = `T2 (hours)`, color = Project, fill = Project)) +
  geom_bar(aes(y = y2 * 60), width = 0.9, alpha = 0.2, position = 'stack', stat = 'identity') +
  geom_point(size = 2) +
  geom_line(aes(group = Project), linewidth = 1) +
  geom_text(aes(label = `T2 (hours)`), size = 5, vjust = -0.5, hjust = 0.5,
            show.legend = F) +
  labs(x = '', y = 'Cumulative time (hours)') +
  scale_x_discrete(expand = c(0, 0), name = '',
                   # sec.axis = sec_axis(~., breaks = num2chr, labels = levels(dfplot$Process))
  )+
  scale_y_continuous(expand = c(0, 0), name = 'Cumulative time (hours)',
                     sec.axis = sec_axis(~./60, name = 'Percent', labels = scales::percent)) +
  # coord_flip() +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15, family = 'Arial Unicode MS'),
        axis.text.x.bottom = element_text(size = 20, family = 'Arial Unicode MS'),
        legend.position = 'top')

dfplot_lgd <- dfplot %>% distinct(x, Process) %>% arrange(x)
p2 <- ggplot(dfplot_lgd, aes(x = Process, color = Process, fill = x)) +
  geom_bar(width = 1) +
  scale_fill_manual(values = rep('white', nrow(dfplot_lgd))) +
  scale_color_manual(values = rep('white', nrow(dfplot_lgd))) +
  labs(x = '', y = 'Cumulative time (hours)') +
  theme_void() +
  theme(text = element_text(size = 15, family = 'Arial Unicode MS'),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = 'left')

p <- ggpubr::ggarrange(p1, p2)

pdf('F6B_timeline_v2.pdf', width = 13, height = 5)
showtext_begin()
print(p)
showtext_end()
graphics.off()

list(Figure6B = dfplot) %>% rio::export('F6B_timeline_v2.xlsx')


# 3.new style -2023.12.29 ------------
mat_procedure <- dfplot %>%
  pivot_wider(id_cols = Project, names_from = 'Process', values_from = `T2 (hours)`) %>% 
  column_to_rownames('Project') %>% 
  .[, levels(dfplot$Process)] %>% 
  add_column('0' = 0, .before = 1)
for(i in 1:nrow(mat_procedure)){
  for(j in 1:ncol(mat_procedure)){
    if(is.na(mat_procedure[i, j])){
      mat_procedure[i, j] <- mat_procedure[i, j-1]
    }
  }
}

dfplot_new <- mat_procedure %>% rownames_to_column('Project') %>% 
  pivot_longer(cols = -Project, names_to = 'Process', values_to = 'y1') %>% 
  as.data.frame()
dfplot_new$Process %<>% factor(levels = c('0', as.character(progress_all)))
dfplot_new$x <- as.numeric(dfplot_new$Process) - 1
dfplot_new[dfplot_new$x != 0, 'x'] <- num2chr[dfplot_new[dfplot_new$x != 0, 'x']]
dfplot_new$x %<>% factor(levels = c('0', num2chr))

p1 <- ggplot(dfplot_new, aes(x = x, y = y1, color = Project, fill = Project)) +
  # geom_bar(aes(y = y2 * 60), width = 0.9, alpha = 0.2, position = 'stack', stat = 'identity') +
  geom_point(size = 2) +
  geom_line(aes(group = Project), linewidth = 1) +
  geom_text(data = dfplot_new %>% filter(y1 != 0),
            aes(label = y1), size = 5, vjust = -0.5, hjust = 0.5,
            show.legend = F) +
  labs(x = 'Procedures', y = 'Cumulative time (hours)') +
  # scale_x_discrete(#expand = c(0, 0),
  #   name = 'Procedures',
  #   # sec.axis = sec_axis(~., breaks = num2chr, labels = levels(dfplot_new$Process))
  # )+
  # scale_y_continuous(#expand = c(0, 0),
  #   name = 'Cumulative time (hours)',
  #   sec.axis = sec_axis(~./60, name = 'Percent', labels = scales::percent)) +
  # coord_flip() +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(text = element_text(size = 15, family = 'Arial Unicode MS'),
        axis.text.x.bottom = element_text(size = 20, family = 'Arial Unicode MS'),
        legend.position = 'top')


p <- ggpubr::ggarrange(p1, p2)
pdf('F6B_timeline_v3.pdf', width = 13, height = 5)
showtext_begin()
print(p)
showtext_end()
graphics.off()

list(Figure6B = dfplot_new) %>% rio::export('F6B_timeline_v3.xlsx')



# save data -------------------------------------------------------
save.image('Figure6.RData')

rm(list = ls())

library(tidyverse)

df_id = data.table::fread("results/step1_get_id/id2species.txt", header = TRUE)

df_pep = data.table::fread("results/step4_seq_related/pep.mw.txt", header = TRUE) %>% 
  dplyr::rename(geneid = `#ID`) %>% 
  merge(df_id, by = "geneid", all.x = TRUE) %>% 
  dplyr::group_by(species) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(n)

# 计数
df_pep %>% 
  dplyr::select(length, species) %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(n = n()) %>% 
  as.data.frame() %>% 
  print()
  
# 常规统计

# 长度anova
df_pep %>% 
  dplyr::select(length, species) -> df.temp

aov(length ~ species, data = df.temp) %>% summary()

df_pep %>% 
  dplyr::select(length, species) %>% 
  dplyr::group_by(species) %>% 
  summarise(min = min(length),
            max = max(length),
            mean = mean(length),
            median = median(length),
            sd = sd(length),
            q1 = quantile(length, 0.25),
            q3 = quantile(length, 0.75)) -> df.summary.length

df.summary.length %>% 
  reshape2::melt(id.vars = 1) %>% 
  ggplot(aes(variable, value, fill = variable)) +
  geom_boxplot() +
  scale_fill_aaas() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")


#####
# 蛋白质量
df_pep %>% 
  dplyr::select(`MW(Da)`, species) -> df.temp

df_pep %>% 
  dplyr::select(`MW(Da)`, species) %>% 
  dplyr::rename(length = `MW(Da)`) %>% 
  dplyr::group_by(species) %>% 
  summarise(min = min(length),
            max = max(length),
            mean = mean(length),
            median = median(length),
            sd = sd(length),
            q1 = quantile(length, 0.25),
            q3 = quantile(length, 0.75)) -> df.summary.mw



# Pi
df_pep %>% 
  dplyr::select(pI, species) -> df.temp

df_pep %>% 
  dplyr::select(pI, species) %>% 
  dplyr::rename(length = pI) %>% 
  dplyr::group_by(species) %>% 
  summarise(min = min(length),
            max = max(length),
            mean = mean(length),
            median = median(length),
            sd = sd(length),
            q1 = quantile(length, 0.25),
            q3 = quantile(length, 0.75)) -> df.summary.pi

library(openxlsx)
names <- list('length' = df.summary.length, 'mw' = df.summary.mw, 'pi' = df.summary.pi)
write.xlsx(names, file = 'results/step4_seq_related/summary.xlsx')



df_pep$species = factor(df_pep$species, levels = unique(df_pep$species))

# 基因蛋白长度
df_pep %>% 
  dplyr::select(length, species) %>% 
  ggplot(aes(x = length, y = species),show.legend = FALSE) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  geom_vline(xintercept = mean(df_pep$length), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(100,1500,100)) +
  labs(x = "Protein length", y = "Rice accession") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"))

export::graph2tif(file = "results/step4_seq_related/蛋白长度统计.tif",
                  width = 10,
                  height = 8,
                  dpi = 500)
# 基因蛋白长度
df_pep %>% 
  dplyr::select(pI, species) %>% 
  ggplot(aes(x = pI, y = species),show.legend = FALSE) +
  geom_boxplot(width = 0.6, outlier.size = 0.8) +
  geom_vline(xintercept = mean(df_pep$pI), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(4,12,1)) +
  labs(x = "Protein isoelectric point", y = "Rice accession") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12, color = "black"))

export::graph2tif(file = "results/step4_seq_related/蛋白分子等电点统计.tif",
                  width = 10,
                  height = 8,
                  dpi = 500)

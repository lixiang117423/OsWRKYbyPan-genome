rm(list = ls())

library(ggsci)
library(ggprism)

all_spe = NULL

df = data.table::fread("results/step2_hmm_search/unique.id.new.wrky.txt", header = FALSE) %>% 
  dplyr::mutate(species = "")

spe = dir("data/geneIDandmRNAid/geneID2mRNAid/")

species = c()

for (i in spe) {
  temp = stringr::str_split(i,"\\.")[[1]][1]
  species = c(species, temp)
}

for (i in 1:nrow(df)) {
  if (stringr::str_sub(df$V1[i],1,3) == "LOC") {
    df$species[i] = "Nipponbare"
  }else{
    for (j in species) {
      if (grepl(j,df$V1[i])) {
        df$species[i] = j
      }
    }
  }
}

colnames(df) = c("geneid","species")

data.table::fwrite(df, file = "results/step1_get_id/id2species.txt",row.names = FALSE, col.names = TRUE)


# 染色体分布偏好性
table(df$species) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mean = mean(Freq)) %>% 
  dplyr::arrange(Freq) -> df.spe

levels = unique(df.spe$Var1)


df_id = data.table::fread("results/step1_get_id/mRNA2geneID.txt", header = TRUE)

df_chr = merge(df, df_id, by.x = "geneid", by.y = "#mRNA_ID", all.x = TRUE) %>% 
  dplyr::group_by(species, chr) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(chr != "Contig109END")
df_chr$chr = factor(df_chr$chr, levels = c(paste0("Chr",1:12),"Contig109END"))

ggplot(df_chr, aes(x = n, y = species, fill = chr)) +
  geom_bar(stat = "identity", width = 1, position="fill") +
  #geom_vline(xintercept = 105, color = "white") +
  scale_fill_igv() +
  scale_x_continuous(expand = c(0,0), label = scales::percent_format(scale = 100)) +
  scale_y_discrete(limit = levels) +
  labs(x = "Percent", y = "Rice Accessions",
       title = "Percent of genes in rice accessions",
       fill = "Chromosome") +
  theme_prism() +
  theme(legend.title = element_text("Chromosome"),
        plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 12)) -> p1
p1
export::graph2tif(file = "results/step4_seq_related/wrky基因数量统计_染色体.tif",
                  width = 8,
                  height = 8,
                  dpi = 500)

# 染色体上基因计数
df_chr %>% 
  dplyr::select(species, chr, n) %>% 
  dplyr::mutate(temp = paste0(species, chr)) %>% 
  dplyr::filter(!duplicated(temp)) %>% 
  dplyr::select(-temp) %>% 
  tidyr::pivot_wider(names_from = chr, values_from = n) %>% 
  tibble::column_to_rownames(var = "species")  %>% 
  pheatmap::pheatmap(cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     angle_col = 45,
                     legend = FALSE,
                     display_numbers = TRUE,
                     number_format = "%.0f",
                     number_color = "black",
                     cellwidth = 20) %>% 
  ggplotify::as.ggplot()

export::graph2tif(file = "results/step4_seq_related/chr.wrky.num.plot.tif",
                  width = 12,
                  height = 6,
                  dpi = 500)

df_chr %>% 
  dplyr::select(species, chr, n) %>% 
  ggplot(aes(x = chr, y = species, fill = n)) +
  geom_tile(aes(width = 1)) +
  geom_text(aes(label = n)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",  
                       midpoint = 8, limit = c(2,24), space = "Lab") +
  labs(x = "Chromosome",
       title = "Number of genes in rice accessions") +
  theme_prism() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12)) -> p2
p2

# 组合图
p2 %>%  aplot::insert_left(p1, width = 0.6)

export::graph2tif(file = "results/step4_seq_related/chr.wrky.num.percent.plot.tif",
                  width = 12,
                  height = 8,
                  dpi = 500)


# 导入序列亚分组
library(ggsci)
library(ggprism)

df_subgroup = data.table::fread("results/step4_seq_related/seq_subgroup_chenxujun.txt", header = TRUE) %>% 
  dplyr::mutate(id = stringr::str_replace(id,">",""))
for (i in 1:nrow(df_subgroup)) {
  df_subgroup$id[i] = stringr::str_split(df_subgroup$id[i]," ")[[1]][1]
}
df_subgroup_plot = merge(df_subgroup, df, by.x = "id", by.y = "V1") %>% 
  dplyr::mutate(temp2 = paste0(species, subgroup)) %>% 
  dplyr::mutate(temp2 = stringr::str_replace_all(temp2," ","")) %>% 
  dplyr::group_by(temp2) %>% 
  dplyr::mutate(n = n(),
                sum = sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(species, subgroup, n, temp2, sum)
  
df_subgroup_plot = df_subgroup_plot[!duplicated(df_subgroup_plot$temp2),] %>% 
  dplyr::arrange(sum)
df_subgroup_plot$species = factor(df_subgroup_plot$species, levels = unique(df_subgroup_plot$species))

table(df$species) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mean = mean(Freq)) %>% 
  dplyr::arrange(Freq) -> df.spe

levels = unique(df.spe$Var1)

ggplot(df_subgroup_plot, aes(x = n, y = species, fill = subgroup)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_vline(xintercept = 105, color = "white") +
  scale_fill_igv() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(limit = levels) +
  labs(x = "Number of Genes", y = "Rice Accessions",
       fill = "Subgroup",
       title = "Number of genes in rice accessions") +
  theme_prism() 

export::graph2tif(file = "results/step4_seq_related/wrky基因数量统计_subgroup.tif",
                  width = 8,
                  height = 8,
                  dpi = 500)


# CH分组
# 导入序列亚分组
library(ggsci)
library(ggprism)

df_subgroup = data.table::fread("results/step4_seq_related/seq_subgroup_CH_type.txt", header = TRUE) %>% 
  dplyr::mutate(id = stringr::str_replace(id,">",""))
for (i in 1:nrow(df_subgroup)) {
  df_subgroup$id[i] = stringr::str_split(df_subgroup$id[i]," ")[[1]][1]
}
df_subgroup_plot = merge(df_subgroup, df, by.x = "id", by.y = "V1") %>% 
  dplyr::mutate(temp2 = paste0(species, subgroup)) %>% 
  dplyr::mutate(temp2 = stringr::str_replace_all(temp2," ","")) %>% 
  dplyr::group_by(temp2) %>% 
  dplyr::mutate(n = n(),
                sum = sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(species, subgroup, n, temp2, sum)

df_subgroup_plot = df_subgroup_plot[!duplicated(df_subgroup_plot$temp2),] %>% 
  dplyr::arrange(sum)
df_subgroup_plot$species = factor(df_subgroup_plot$species, levels = unique(df_subgroup_plot$species))

table(df$species) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mean = mean(Freq)) %>% 
  dplyr::arrange(Freq) -> df.spe

levels = unique(df.spe$Var1)

ggplot(df_subgroup_plot, aes(x = n, y = species, fill = subgroup)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_vline(xintercept = 105, color = "white") +
  scale_fill_igv() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(limit = levels) +
  labs(x = "Number of Genes", y = "Rice Accessions",
       fill = "Subgroup",
       title = "Number of genes in rice accessions") +
  theme_prism() 

export::graph2tif(file = "results/step4_seq_related/wrky基因数量统计_subgroup_CH.tif",
                  width = 8,
                  height = 8,
                  dpi = 500)



# 统计每个species的数量
library(ggprism)
library(ggplot2)

table(df$species) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mean = mean(Freq)) %>% 
  dplyr::arrange(Freq) -> df.spe
  
df.spe$Var1 = factor(df.spe$Var1, levels = unique(df.spe$Var1))

all_spe = df.spe %>% 
  dplyr::rename(pfam = Freq) %>% 
  dplyr::select(1:2)

ggplot(df.spe) +
  geom_bar(aes(y = Var1, x = Freq), stat = "identity", width = 0.8) +
  geom_text(aes(x = Freq + 2, y = Var1, label = Freq)) +
  geom_vline(xintercept = 97, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 110, color = "white") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Number of Genes", y = "Rice Accessions", title = "Number of genes in rice accessions") +
  theme_prism()

export::graph2tif(file = "results/step2_hmm_search/wrky基因数量统计_pfam.tif",
                  width = 8,
                  height = 8,
                  dpi = 500)

# planttdb indica
rm(list = ls())

df = data.table::fread("results/step3_planttfdb_rebuild_and_research/planttfdb.indica.unique.id.txt", header = FALSE) %>% 
  dplyr::mutate(species = "")

spe = dir("data/geneIDandmRNAid/geneID2mRNAid/")

species = c()

for (i in spe) {
  temp = stringr::str_split(i,"\\.")[[1]][1]
  species = c(species, temp)
}

for (i in 1:nrow(df)) {
  if (stringr::str_sub(df$V1[i],1,3) == "LOC") {
    df$species[i] = "Nipponbare"
  }else{
    for (j in species) {
      if (grepl(j,df$V1[i])) {
        df$species[i] = j
      }
    }
  }
}


# 统计每个species的数量

table(df$species) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mean = mean(Freq)) %>% 
  dplyr::arrange(Freq) -> df.spe

df.spe$Var1 = factor(df.spe$Var1, levels = unique(df.spe$Var1))


all_spe = merge(all_spe, df.spe[,c(1,2)], by = "Var1")
colnames(all_spe)[3] = "PlantTFDB_indica"

ggplot(df.spe) +
  geom_bar(aes(y = Var1, x = Freq), stat = "identity", width = 0.8) +
  geom_text(aes(x = Freq + 3, y = Var1, label = Freq)) +
  geom_vline(xintercept = 97, linetype = "dashed") +
  geom_vline(xintercept = 110, color = "white") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Number of Genes", y = "Rice Accessions", title = "Number of genes in rice accessions (PlantTFDB indica)") +
  theme_prism()

export::graph2tif(file = "results/step2_hmm_search/wrky基因数量统计_PlantTFDB indica.tif",
                  width = 8,
                  height = 8,
                  dpi = 500)




# planttdb japonica
rm(list = ls())

df = data.table::fread("results/step3_planttfdb_rebuild_and_research/planttfdb.japonica.unique.id.txt", header = FALSE) %>% 
  dplyr::mutate(species = "")

spe = dir("data/geneIDandmRNAid/geneID2mRNAid/")

species = c()

for (i in spe) {
  temp = stringr::str_split(i,"\\.")[[1]][1]
  species = c(species, temp)
}

for (i in 1:nrow(df)) {
  if (stringr::str_sub(df$V1[i],1,3) == "LOC") {
    df$species[i] = "Nipponbare"
  }else{
    for (j in species) {
      if (grepl(j,df$V1[i])) {
        df$species[i] = j
      }
    }
  }
}


# 统计每个species的数量

table(df$species) %>% 
  as.data.frame() %>% 
  dplyr::mutate(mean = mean(Freq)) %>% 
  dplyr::arrange(Freq) -> df.spe

df.spe$Var1 = factor(df.spe$Var1, levels = unique(df.spe$Var1))

all_spe = merge(all_spe, df.spe[,c(1,2)], by = "Var1")
colnames(all_spe)[4] = "PlantTFDB_japonica"

ggplot(df.spe) +
  geom_bar(aes(y = Var1, x = Freq), stat = "identity", width = 0.8) +
  geom_text(aes(x = Freq + 3, y = Var1, label = Freq)) +
  geom_vline(xintercept = 97, linetype = "dashed") +
  geom_vline(xintercept = 110, color = "white") +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Number of Genes", y = "Rice Accessions", title = "Number of genes in rice accessions (PlantTFDB japonica)") +
  theme_prism()

export::graph2tif(file = "results/step2_hmm_search/wrky基因数量统计_PlantTFDB japonica.tif",
                  width = 8,
                  height = 8,
                  dpi = 500)


























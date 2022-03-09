rm(list = ls())

library(ggtree)
library(ggtreeExtra)
library(ggsci)
library(ggnewscale)

df_tree = ggtree::read.tree("results/step4_seq_related/new.id.ch.type.aligned.tree.nwk")


# 导入序列亚群分组

# CH类型分组
df_seq_group = data.table::fread("results/step4_seq_related/seq_subgroup_CH_type.txt", header = TRUE) %>% 
  dplyr::mutate(ID = id, group = subgroup) %>% 
  dplyr::mutate(ID = stringr::str_replace(id,">","")) %>% 
  dplyr::select(ID, group)

# 陈旭君等分组
df_seq_subgroup = data.table::fread("results/step4_seq_related/seq_subgroup_chenxujun.txt", header = TRUE) %>% 
  dplyr::select(id, subgroup) %>% 
  dplyr::rename(ID = id) %>% 
  dplyr::mutate(ID = stringr::str_sub(ID,2, nchar(ID)))

# 7肽分组
df_7_tai = data.table::fread("results/step4_seq_related/7taigroup.txt", header = TRUE) %>% 
  dplyr::rename(ID = id,
                WRKY = first_char) %>% 
  dplyr::mutate(ID = stringr::str_replace(ID,">",""))
table(df_7_tai$WRKY) %>% as.data.frame() %>% dplyr::summarise(sum = sum(Freq))

# NCBI-CDD pvalue
df_ncbi = data.table::fread("results/step4_seq_related/ncbi_cdd_search.txt", header = TRUE)

for (i in 1:nrow(df_ncbi)) {
  df_ncbi$Query[i] = stringr::str_split(df_ncbi$Query[i],">")[[1]][2]
}

df_ncbi_selected = df_ncbi %>% 
  dplyr::filter(Query %in% df_seq_group$ID) %>% 
  #dplyr::filter(`Short name` == "WRKY") %>% 
  dplyr::filter(!duplicated(Query)) %>% 
  dplyr::rename(ID = Query) %>% 
  merge(df_seq_group, by = "ID")

# 处理第一个结构域序列
df_domain_1 = pac4xiang::fasta2df("results/step2_hmm_search/first.domain.fa") %>% 
  dplyr::mutate(id = stringr::str_replace(id, ">", ""))

for (i in 1:nrow(df_domain_1)) {
  df_domain_1$id[i] = stringr::str_split(df_domain_1$id[i]," ")[[1]][1]
}

df_domain_1 %>% 
  dplyr::filter(id %in% df_seq_group$ID) %>% 
  pac4xiang::df2fasta() %>% 
  data.table::fwrite(file = "results/step4_seq_related/first_domain.fa",
                     col.names = FALSE, row.names = FALSE, quote = FALSE)
  



# 绘制进化树
p = ggtree(df_tree, layout = "fan", open.angle = 0, branch.length = "none")
p

p1 = p %<+% df_seq_subgroup +
  geom_tippoint(aes(color = subgroup)) +
  scale_color_igv(name = "Subgroup")
p1


p2 = p1 %<+%  df_seq_group + 
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = group),
             offset = 0.04, size = 0.02, width = 1.5
             ) +
  scale_fill_aaas(name = " Zinc Finger Type")

p2

# 20220308 add
p3 = p2 + new_scale_fill() +
  geom_fruit(data = df_7_tai,
             geom = geom_tile,
             mapping = aes(y = ID,fill = WRKY),
             offset = 0.05, size = 0.02, width = 1.5
  ) +
  scale_fill_lancet(name = "WRKY Type")

p3

# 区分亚洲和非洲
df_area = df_seq_group %>% 
  dplyr::select(ID) %>% 
  dplyr::mutate(area = stringr::str_extract(ID,"CG14")) %>% 
  dplyr::mutate(area = case_when(area == "CG14" ~ "Africa",
                                 is.na(area) ~ "Asia"))

p4 = p3 + new_scale_fill() +
  geom_fruit(data = df_area,
             geom = geom_tile,
             mapping = aes(y = ID,fill = area),
             offset = 0.05, size = 0.02, width = 1.5
  ) +
  scale_fill_nejm(name = "Region")

p4

export::graph2tif(file = "results/step4_seq_related/phylo.tree.tif",
                  width = 10, height = 8, dpi = 500)

################################################################################
# 匹配序列
p = ggtree(df_tree, layout = "fan", open.angle = 0)
p

p1 = p %<+% df_seq_subgroup +
  geom_tippoint(aes(color = subgroup)) +
  scale_color_igv(name = "Subgroup")
p1


p2 = p1 %<+%  df_seq_group + 
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = group),
             offset = 0.04, size = 0.02, width = 0.3
  ) +
  scale_fill_aaas(name = "CH Type")

p2

p3 = p2 +  layout_rectangular() +
  theme(legend.position = c(0.15,0.6),
        legend.text=element_text(size=13),
        legend.title = element_text(size = 15))
p3

export::graph2tif(file = "results/step4_seq_related/phylo.tree.水平.tif",
                  width = 10, height = 8, dpi = 500)






rm(list = ls())

library(ggplot2)
library(ggsci)

df1 = data.table::fread("results/LOC_Os02g08440表达量-南方科技大学网站/LOC_OS02G08440_biotic_expression_levels.csv", header = TRUE)

ggplot(df1, aes(x = biotic, y = LOC_OS02G08440, fill = biotic)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  coord_flip() +
  scale_fill_igv() +
  scale_y_continuous(breaks = seq(0,300,50)) +
  labs(title = "Expression Level under Biotic Stress") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12)) -> p1
p1

# 非生物胁迫
df2 = data.table::fread("results/LOC_Os02g08440表达量-南方科技大学网站/LOC_OS02G08440_abiotic_expression_levels.csv" ,header = TRUE)

ggplot(df2, aes(x = abiotic, y = LOC_OS02G08440, fill = abiotic)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  coord_flip() +
  scale_fill_igv() +
  scale_y_continuous(breaks = seq(0,400,50)) +
  labs(title = "Expression Level under Abiotic Stress") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12)) -> p2
p2


# 不同部位
df3 = data.table::fread("results/LOC_Os02g08440表达量-南方科技大学网站/LOC_OS02G08440_tissue_expression.csv", header = TRUE)

ggplot(df3, aes(x = `Tissue/Satge`, y = LOC_OS02G08440, fill = `Tissue/Satge`)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  coord_flip() +
  scale_fill_igv() +
  scale_y_continuous(breaks = seq(0,400,50)) +
  labs(title = "Expression Level at Different Tissue") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12)) -> p3
p3


# 不同部位
df4 = data.table::fread("results/LOC_Os02g08440表达量-南方科技大学网站/LOC_OS02G08440_tissue_expression (1).csv", header = TRUE)
ggplot(df4, aes(x = `Tissue/Satge`, y = LOC_OS02G08440, fill = `Tissue/Satge`)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  coord_flip() +
  scale_fill_igv() +
  scale_y_continuous(breaks = seq(0,300,50)) +
  labs(title = "Expression Level at Different Development Stage") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(colour = "black", size = 12)) -> p4
p4

# combian
p1 + p2 + p3 + p4 + patchwork::plot_layout(ncol = 2)

export::graph2tif(file = "results/各种表达量.tif",
                  width = 20, height = 15, dpi = 500)
export::graph2pdf(file = "results/LOC_Os02g08440表达量-南方科技大学网站/pre.pdf",
                  width = 20, height = 15)

cowplot::plot_grid(p1, p2, p3, p4,
                   labels = LETTERS[1:4],
                   label_size = 25,
                   ncol = 2)
export::graph2tif(file = "results/各种表达量.tif",
                  width = 20, height = 15, dpi = 500)
export::graph2pdf(file = "results/LOC_Os02g08440表达量-南方科技大学网站/pre.pdf",
                  width = 20, height = 15)


# 全部表达量
df_id = data.table::fread("rna-seq/stress/results/DEGs.wrky.txt", header = TRUE) %>% 
  dplyr::select(row) %>% 
  dplyr::distinct_all() %>% 
  dplyr::mutate(row = stringr::str_replace(row, "Os","OS"),
                row = stringr::str_replace(row, "g","G"))

group = data.table::fread("rna-seq/group.csv", header = TRUE) %>% 
  dplyr::filter(Cultivar == "Nipponbare") %>% 
  dplyr::select(Sample) %>% 
  dplyr::distinct_all()

df = data.table::fread("rna-seq/rice_fpkmAllData.csv", header = TRUE) %>% 
  dplyr::filter(Sample %in% df_id$row) %>% 
  tibble::column_to_rownames(var = "Sample") %>% 
  dplyr::select(group$Sample) %>% 
  pheatmap:::scale_rows()

df[df>3.17] = 3.17


pheatmap::pheatmap(df,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   color=colorRampPalette(c("navy", "white", "red"))(50)) ->p

# 提取聚类结果
index <- cutree(p$tree_col,k=2) %>% 
  as.data.frame()
colnames(index) = "tree"
index$sample = rownames(index)
index.1 = index %>% 
  dplyr::filter(tree == 2)


export::graph2tif(file = "rna-seq/pprd.wrky.all.tif",
                  width = 10, height = 10, dpi = 500)


# cor 
df = data.table::fread("rna-seq/rice_fpkmAllData.csv", header = TRUE) %>% 
  dplyr::filter(Sample %in% df_id$row) %>% 
  tibble::column_to_rownames(var = "Sample") %>% 
  dplyr::select(group$Sample) %>% 
  t() %>% 
  as.data.frame()

WGCNA::corAndPvalue(df)$cor %>% 
  as.matrix() %>% 
  corrplot::corrplot()

library(linkET)
library(ggtext)

correlate(df) %>% 
  as_md_tbl() %>% 
  qcorrplot(type = "lower") +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  theme(axis.text = element_blank(),
        axis.title = element_blank())

# 
df.cor = WGCNA::corAndPvalue(df)$cor %>% 
  as.data.frame() %>% 
  dplyr::mutate(temp = 1) %>% 
  reshape2::melt(id.vars = "temp") 



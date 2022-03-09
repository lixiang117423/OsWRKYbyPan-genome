rm(list = ls())

df_id = data.table::fread("results/step2_hmm_search/unique.id.new.wrky.txt",header = FALSE) %>% 
  dplyr::mutate(temp = stringr::str_sub(V1, 3, 6)) %>% 
  dplyr::filter(temp == "IR64") %>% 
  dplyr::mutate(gene = stringr::str_sub(V1, 1, 20)) %>% 
  dplyr::select(gene)

df_smaple = data.table::fread("rna-seq/warmnighttemp/data/sample.info.txt", header = TRUE) %>% 
  dplyr::filter(group1 %in% c("Control"))

df_fpkm = data.table::fread("rna-seq/warmnighttemp/FPKM.csv", header = TRUE) %>% 
  dplyr::filter(V1 %in% df_id$gene) %>% 
  dplyr::mutate(V1 = stringr::str_sub(V1, 1, 17)) %>% 
  tibble::column_to_rownames(var = "V1") %>% 
  dplyr::select(df_smaple$id)


df_ann_col = df_smaple %>% 
  dplyr::select(1,4) %>% 
  tibble::column_to_rownames(var = "id") %>% 
  dplyr::rename(Time = group3) %>% 
  dplyr::mutate(Time = factor(Time, levels = c(0,3.5,7, 10.5, 12, 14, 17.5, 23)))


df_fpkm %>% 
  pheatmap:::scale_rows() %>% 
  t() %>% 
  pheatmap::pheatmap(annotation_row = df_ann_col,
                     show_colnames = FALSE,
                     show_rownames = FALSE,
                     cellwidth = 5,
                     cellheight = 5) -> p.temp
p = ggplotify::as.ggplot(p.temp)

ggsave(p, filename = "results/step7_rna-seq/夜间高温CK时序.png",
       width = 12, height = 6, dpi = 500)


export::graph2ppt(p, file = "results/step7_rna-seq/夜间高温CK时序.pptx",
                  width = 12, height = 6, center = TRUE)






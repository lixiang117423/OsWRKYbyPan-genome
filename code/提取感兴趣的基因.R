rm(list = ls())

library(ggplotify)

df_id = data.table::fread("results/step1_get_id/id2species.txt", header = TRUE) %>% 
  dplyr::filter(species == "IR64") %>% 
  dplyr::mutate(geneid = stringr::str_replace(geneid,".T01",""))

df_sample_info = data.table::fread("rna-seq/data/sample.info.txt", header = TRUE)

df_tpm = data.table::fread("rna-seq/TPM.csv", header = TRUE) %>% 
  dplyr::rename(geneid = V1) %>% 
  dplyr::filter(geneid %in% df_id$geneid) %>% 
  as.data.frame()

rownames(df_tpm) = df_tpm$geneid
df_tpm = df_tpm[,-1] 


df_ann = df_sample_info %>% 
  dplyr::mutate(group3 = as.character(group3)) %>% 
  as.data.frame()
rownames(df_ann) = df_ann$id
df_ann = df_ann[,-1]

df_tpm %>% 
  pheatmap:::scale_rows() %>% 
  pheatmap::pheatmap(show_rownames = FALSE,
                     show_colnames = FALSE,
                     cluster_cols = FALSE,
                     annotation_col = df_ann) ->p 
p
p1 = ggplotify::as.ggplot(p)
p1
export::graph2tif(p1, file = "rna-seq/wrky.heatmap.tif",width = 10, height = 8, dpi = 500)

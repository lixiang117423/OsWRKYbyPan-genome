###################################
# 处理ID
###################################

# indica
df_indica <- pac4xiang::fasta2df("results/step3_planttfdb_rebuild_and_research/wrky_planttfdb_indica.fas")

for (i in 1:nrow(df_indica)) {
  df_indica$id[i] <- stringr::str_split(df_indica$id[i], "\\|")[[1]][1] %>%
    stringr::str_replace(">", "")
}

pac4xiang::df2fasta(df_indica) %>%
  data.table::fwrite(
    file = "results/step3_planttfdb_rebuild_and_research/wrky_planttfdb_indica.new.fas",
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )


# japonica
df_japonica <- pac4xiang::fasta2df("results/step3_planttfdb_rebuild_and_research/wrky_planttfdb_japonica.fas")

for (i in 1:nrow(df_japonica)) {
  df_japonica$id[i] <- stringr::str_split(df_japonica$id[i], "\\|")[[1]][1] %>%
    stringr::str_replace(">", "")
}

pac4xiang::df2fasta(df_japonica) %>%
  data.table::fwrite(
    file = "results/step3_planttfdb_rebuild_and_research/wrky_planttfdb_japonica.new.fas",
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )


# 比较三种方法得到的基因ID的异同

library(ggvenn)

geneid = data.table::fread("results/step1_get_id/mRNA2geneID.txt", header = TRUE) 

df1 <- data.table::fread("results/step2_hmm_search/unique.id.txt", header = FALSE) %>% 
  merge(geneid[,c(1,2)], by.x = "V1", by.y = "#mRNA_ID", all.x = TRUE)

df2 <- data.table::fread("results/step3_planttfdb_rebuild_and_research/planttfdb.indica.unique.id.txt", header = FALSE)  %>% 
  merge(geneid[,c(1,2)], by.x = "V1", by.y = "#mRNA_ID", all.x = TRUE)

df3 <- data.table::fread("results/step3_planttfdb_rebuild_and_research/planttfdb.japonica.unique.id.txt", header = FALSE) %>% 
  merge(geneid[,c(1,2)], by.x = "V1", by.y = "#mRNA_ID", all.x = TRUE)


li <- list(
  PlantTFDB_indica = df2$gene_ID,
  PlantTFDB_japonica = df3$gene_ID,
  Pfam = df1$gene_ID
)

ggvenn(li, 
       show_percentage = FALSE,
       fill_color = c("#ff5900","#0d0dbb","#00b76d"),  
       fill_alpha = 0.8,
       set_name_color = "black",
       set_name_size = 6,
       text_color = "white",
       text_size = 6
       )
export::graph2tif(file = "results/step3_planttfdb_rebuild_and_research/三种方法比较的韦恩图.tif",
                  width = 6,
                  height = 6,
                  dpi = 500)












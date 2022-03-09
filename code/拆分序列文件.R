rm(list = ls())

df_seq = pac4xiang::fasta2df("results/step4_seq_related/extracted.pep.fa") %>% 
  dplyr::mutate(id = stringr::str_replace(id, ">","")) 

for (i in 1:nrow(df_seq)) {
  df_seq$id[i] = stringr::str_split(df_seq$id[i]," ")[[1]][1]
}

df_seq %>% 
  pac4xiang::df2fasta() %>% 
  data.table::fwrite(file = "results/step4_seq_related/extracted.pep.new.id.fa", col.names = FALSE, row.names = FALSE)

for (i in seq(1,3500,500)) {
  if (i == 3001) {
    df_sub = df_seq[3001:3257,] %>% 
      pac4xiang::df2fasta()
    filename = "results/step4_seq_related/seq_500_in_one_file/3001-3257.fa"
    data.table::fwrite(df_sub, file = filename, col.names = FALSE, row.names = FALSE)
  }else{
    df_sub = df_seq[i:(i+499),] %>% 
      pac4xiang::df2fasta()
    filename = paste0("results/step4_seq_related/seq_500_in_one_file/",i,"-",i+499,".fa")
    data.table::fwrite(df_sub, file = filename, col.names = FALSE, row.names = FALSE)
  }
}

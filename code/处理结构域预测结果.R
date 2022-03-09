rm(list = ls())


df_all_id = data.table::fread("results/step2_hmm_search/unique.id.txt",header = FALSE)

# NCBI
df_ncbi = data.table::fread("results/step4_seq_related/seq_500_in_one_file/ncbi_cdd_search.txt", header = TRUE) %>% 
  dplyr::filter(`Hit type` == "specific") %>% 
  dplyr::group_by(Query) %>% 
  dplyr::mutate(n = n()) %>% 
  dplyr::ungroup()

for (i in 1:nrow(df_ncbi)) {
  df_ncbi$Query[i] = stringr::str_split(df_ncbi$Query[i],">")[[1]][2]
}

df_ncbi_1 = df_ncbi %>% 
  dplyr::filter(`Short name` == "WRKY")

df_ncbi_1 = df_ncbi_1[!duplicated(df_ncbi_1$Query),] %>% 
  dplyr::select(1) %>% 
  dplyr::mutate(isin = 1)
colnames(df_ncbi_1)[1] = "V1"

df_ncbi_2 = df_ncbi %>% 
  dplyr::filter(`Short name` != "WRKY")

df_ncbi_not = setdiff(df_all_id$V1, df_ncbi_1$Query) %>% 
  as_data_frame()
colnames(df_ncbi_not) = "V1"
df_ncbi_not = df_ncbi_not %>% 
  dplyr::mutate(isin = 0)

df_ncbi_final = rbind(df_ncbi_1, df_ncbi_not)

# SMART
df_smart_not = data.table::fread("results/step4_seq_related/SMART_not_matched.txt", header = FALSE) %>% 
  dplyr::mutate(isin = 0)

df_smart_is = setdiff(df_all_id$V1, df_smart_not$V1) %>% 
  as.data.frame()
colnames(df_smart_is) = "V1"

df_smart_is = df_smart_is %>% 
  dplyr::mutate(isin = 1)

df_smart = rbind(df_smart_is, df_smart_not)

# Pfam


# merge
df_final = merge(df_ncbi_final, df_smart, by = "V1") %>% 
  merge(df_pfam, by = "V1")
colnames(df_final) = c("id","NCBI-CCD","SMART","Pfam")













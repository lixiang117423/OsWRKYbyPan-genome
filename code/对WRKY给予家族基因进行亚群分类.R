rm(list = ls())

library(tidyverse)

df <- pac4xiang::fasta2df("results/step4_seq_related/new.id.ch.type.fa")

# Paul J. Rushton
df_new <- df %>%
  dplyr::mutate(
    #first_char = stringr::str_extract(seq, "\\w{1}(WRKYGQK)"),
    first_char = stringr::str_extract(seq, "(WRKY)\\w{3}"),
    sec_char = stringr::str_extract(seq, "(RSYYKC)")
  ) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(RSYYKC)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(RSYYRC)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(RGYYKC)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(RGYYRC)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(RAYYRC)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(RSYYRC)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(CPVKKKV)"),
    sec_char
  )) %>%
  dplyr::mutate(sec_char = ifelse(is.na(sec_char),
    stringr::str_extract(seq, "(CPVRKQV)"),
    sec_char
  )) %>%
  dplyr::mutate(
    temp = paste0(first_char, sec_char),
    subgroup = NA
  ) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RSYYKC)")),
    subgroup,
    "Ⅰ NT"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RSYYKC)")),
    subgroup,
    "I CT"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RSYYRC)")),
    subgroup,
    "Ⅱ c"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RGYYKC)")),
    subgroup,
    "Ⅱ d"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RGYYRC)")),
    subgroup,
    "Ⅱ e"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RAYYRC)")),
    subgroup,
    "Ⅲ a"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(RSYYRC)")),
    subgroup,
    "Ⅲ b"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(CPVKKKV)")),
    subgroup,
    "Ⅱ a"
  )) %>%
  dplyr::mutate(subgroup = ifelse(is.na(stringr::str_extract(temp, "(WRKYGQK)(CPVRKQV)")),
    subgroup,
    "Ⅱ b"
  )) %>% 
  dplyr::mutate(subgroup = ifelse(is.na(subgroup),"other",subgroup))

# 输出7肽分组信息
df_new %>% 
  dplyr::select(id, first_char) %>% 
  data.table::fwrite(file = "results/step4_seq_related/7taigroup.txt", row.names = FALSE)


data.table::fwrite(df_new, file = "results/step4_seq_related/seq_subgroup_chenxujun.txt", row.names = FALSE)

# Paul J. Rushton
df_new <- df %>%
  dplyr::mutate(
    #first_char = stringr::str_extract(seq, "\\w{1}(WRKYGQK)"),
    first_char = stringr::str_extract(seq, "(WRKY)"),
    sec_char = stringr::str_extract(seq, "(C)\\w{4,7}(C)"),
    third_char = stringr::str_extract(seq, "(H)\\w{1}(H)")
  ) %>% 
  dplyr::mutate(third_char = ifelse(is.na(third_char),
                                    stringr::str_extract(seq, "(H)\\w{1}(C)"),
                                    third_char
                                    )) %>% 
  dplyr::mutate(lab1 = paste0(stringr::str_sub(sec_char,1,1),
                              stringr::str_sub(sec_char,nchar(sec_char),nchar(sec_char))),
                lab2 = paste0(stringr::str_sub(third_char,1,1),
                              stringr::str_sub(third_char,nchar(third_char),nchar(third_char)))
                ) %>% 
  dplyr::mutate(temp = paste0(lab1, lab2)) %>% 
  dplyr::mutate(subgroup = case_when(temp == "CCHH" ~ "C2H2",
                                  temp == "CCHC" ~ "C2HC",
                                  temp == "CCNANA" ~ "C2-",
                                  temp == "NANAHC" ~ "-HC",
                                  temp == "NANAHH" ~ "-H2",
                                  TRUE ~ "-"))

data.table::fwrite(df_new, file = "results/step4_seq_related/seq_subgroup_CH_type.txt", row.names = FALSE)



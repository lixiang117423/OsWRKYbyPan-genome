rm(list = ls())

library(tidyverse)

data.table::fread("rna-seq/blast/data/genome/MSU.IGDBv1.Allset.gff", header = FALSE) %>% 
  dplyr::mutate(V9= paste0(V12,V13,";")) %>% 
  dplyr::select(1:9) %>% 
  data.table::fwrite(file = "rna-seq/data/genome/ir64.new.1.gtf", 
                     row.names = FALSE, col.names = FALSE, sep = "\t")


data.table::fread("rna-seq/blast/data/genome/MSU.IGDBv1.Allset.gff", header = FALSE) %>% 
  dplyr::mutate(V1= paste0(V1,"_Nip")) %>% 
  #dplyr::select(1:8) %>% 
  data.table::fwrite(file = "rna-seq/blast/data/genome/Nip.gff", 
                     row.names = FALSE, col.names = FALSE, sep = "\t")

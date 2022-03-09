rm(list = ls())

df = data.table::fread("rna-seq/mapping.rate.stat", header = FALSE)

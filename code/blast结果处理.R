rm(list = ls())

library(tidyverse)

df = data.table::fread("results/step6_blast_and_function/blast.res.txt", header = FALSE) %>% 
  dplyr::filter(V3 == 100)


for (i in unique(df$V2)) {
  print(i)
}


library(ggmsa)

ggmsa("results/step6_blast_and_function/gene.pep.aligned.fa",color = "Chemistry_AA") +
  geom_seqlogo() + 
  #geom_msaBar()
  theme

export::graph2png(file = "results/step6_blast_and_function/比对结果可视化.png",
                  width = 60, height = 20, dpi = 500)

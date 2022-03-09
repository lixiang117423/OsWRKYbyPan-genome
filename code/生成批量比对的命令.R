library(tidyverse)

dir_fastq = dir("/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/data/fastq/")

df = NULL

for (i in dir_fastq) {
  if (i == "code") {
    next
  }else{
    j = stringr::str_split(i, "\\.")[[1]][1]
    comm = paste0("hisat2 --new-summary -p 70 -x data/genome/index/ir64.index -U data/fastq/",
                  i," -S mapping/",j,".sam"," R > mapping.log/", j,".mapping.log  2>&1") %>% 
      as.data.frame()
    df = rbind(df, comm)
  }
}


data.table::fwrite(df,
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/run.mapping.sh", 
                   col.names = FALSE, row.names = FALSE)


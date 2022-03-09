rm(list = ls())

library(tidyverse)

df = NULL

dir_fastq = dir("/home/lixiang/project/acuce/rna-seq-15_69_87/data/raw_data/")

sample = c()

for (i in dir_fastq) {
  sample_temp = stringr::str_sub(i,1,8)
  sample = c(sample,sample_temp)
}

sample = unique(sample) %>% 
  setdiff("md5.txt")

for (i in sample) {
  
  if (stringr::str_sub(i, 1, 3) == "A87") {
    comm = paste0("hisat2 --new-summary -p 70 -x", 
                  " /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/genome/index/nip.index",
                  " --min-intronlen 20 --max-intronlen 10000 --dta ",
                  " -1 /home/lixiang/project/acuce/rna-seq-15_69_87/data/raw_data/",i,"_1.fq.gz",
                  " -2 /home/lixiang/project/acuce/rna-seq-15_69_87/data/raw_data/",i,"_2.fq.gz",
                  " -S /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/plant-plant/mapping/",i,".sam") %>% 
      as.data.frame()
    df = rbind(df, comm)
    
    comm1 = paste0("samtools sort -@ 70 -o /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/plant-plant/bam/",i,".sorted.bam", 
                   " /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/plant-plant/mapping/",i,".sam") %>% 
      as.data.frame()
    df = rbind(df, comm1)
    
    comm2 = paste0("samtools index /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/plant-plant/bam/",i,".sorted.bam") %>% 
      as.data.frame()
    df = rbind(df, comm2)
    
    
    comm2_1 = paste0("stringtie ","/home/lixiang/project/acuce/lincRNA_test/bam/",i,".sorted.bam -G ",
                     "/home/publicdata/refgenome/rice/acuce/annotation/YLG.gff3  -j 3 -o ",
                     "/home/lixiang/project/acuce/lincRNA_test/stringtie/",i,".gtf") %>% 
      as.data.frame()
    
    #df = rbind(df, comm2_1)
    
    if (i == "A87_87_3343") {
      comm = paste0("stringtie --merge -F 0.5 -T 0.5 -G ",
                    "/home/publicdata/refgenome/rice/acuce/annotation/YLG.gff3 " ,
                    paste0("/home/lixiang/project/lincRNA/stringtie/",sample,".gtf"),
                    " > Osat.merged.gtf")
    }
  }
}


data.table::fwrite(df,
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/plant-plant/mapping.sh", 
                   col.names = FALSE, row.names = FALSE, quote = FALSE)
 












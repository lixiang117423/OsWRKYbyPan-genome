library(tidyverse)

df = NULL

dir_fastq = dir("/home/lixiang/project/acuce/rna-seq-15_69_87/data/raw_data/")

sample_id = c()
for (i in dir_fastq) {
  if (i == "md5.txt") {
    next
  }else{
    j = stringr::str_sub(i, 1,8)
    sample_id = c(sample_id, j)
  }
}


sample_id = unique(sample_id)

for (i in sample_id) {
  if (i == "md5.txt") {
    next
  }else{
    comm = paste0("hisat2 --new-summary -p 75 -x /opt/publicdata/refgenome/rice/acuce/index/hisat2/acuce.index",
                  " -1 /home/lixiang/project/acuce/rna-seq-15_69_87/data/raw_data/",i,"_1.fq.gz",
                  " -2 /home/lixiang/project/acuce/rna-seq-15_69_87/data/raw_data/",i,"_2.fq.gz",
                  " -S /home/lixiang/project/acuce/rna-seq-15_69_87/mapping/",i,".sam") %>% 
      as.data.frame()
    df = rbind(df, comm)
  
    comm1 = paste0("samtools sort -@ 70 -o /home/lixiang/project/acuce/rna-seq-15_69_87/bam/",i,".sorted.bam", 
                   " /home/lixiang/project/acuce/rna-seq-15_69_87/mapping/",i,".sam") %>% 
      as.data.frame()
    df = rbind(df, comm1)
    
    comm2 = paste0("samtools index /home/lixiang/project/acuce/rna-seq-15_69_87/bam/",i,".sorted.bam") %>% 
      as.data.frame()
    df = rbind(df, comm2)
    
    comm3 = paste0("htseq-count -f sam -r name -i gene_id -s reverse -t exon -m intersection-nonempty",
                   ' /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/mapping/' , i, #, '.sorted.bam' , 
                   ' /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/data/genome/ir64.gtf',
                   ' > ',  '/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/counts/',  j, '.count.txt') %>% 
      as.data.frame()
    
    comm4 = paste0("featureCounts -p -t exon -g gene_id -a /opt/publicdata/refgenome/rice/acuce/annotation/acuce.gtf -o ",
                   "/home/lixiang/project/acuce/rna-seq-15_69_87/counts/",i,".counts.featurecounts.txt ",
                   "/home/lixiang/project/acuce/rna-seq-15_69_87/bam/",i,".sorted.bam")  %>% 
      as.data.frame()
    df = rbind(df, comm4)
    
    comm5 = paste0("cut -f 1,7 /home/lixiang/project/acuce/rna-seq-15_69_87/counts/",i,".counts.featurecounts.txt | grep -v '^#' > ", 
                   "/home/lixiang/project/acuce/rna-seq-15_69_87/counts/",i,".counts.txt") %>% 
      as.data.frame()
    df = rbind(df, comm5)
  }
}


data.table::fwrite(df,
                   "/home/lixiang/project/acuce/rna-seq-15_69_87/mapping.sam2bam.and.featurecounts.sh", 
                   col.names = FALSE, row.names = FALSE, quote = FALSE)













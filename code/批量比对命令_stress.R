library(tidyverse)

df = NULL

dir_fastq = dir("/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/rawdata/fastq/")

sample = data.frame()

for (i in dir_fastq) {
  fa_temp = stringr::str_split(i, "\\.")
  len = length(stringr::str_split(fa_temp[[1]][1],"_")[[1]])
  if (len == 1) {
    df_temp = data.frame(sra = stringr::str_split(fa_temp[[1]][1],"_")[[1]][1],
                         group = "single")
    
  }else{
    df_temp = data.frame(sra = stringr::str_split(fa_temp[[1]][1],"_")[[1]][1],
                         group = "pair")
    
  }
  
  sample = rbind(sample, df_temp)
}

sample = sample %>% 
  dplyr::filter(!duplicated(sra))


for (i in 1:nrow(sample)) {
  
  m = sample$sra[i]
  
  if (sample$group[i] == "single") {
    comm = paste0("hisat2 --new-summary -p 75 -x /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/genome/index/nip.index",
                  " -U /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/rawdata/fastq/",m,".fastq",
                  " -S /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/mapping/",m,".sam") %>% 
      as.data.frame()
    df = rbind(df, comm)
  }else{
    comm = paste0("hisat2 --new-summary -p 75 -x /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/genome/index/nip.index",
                  " -1 /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/rawdata/fastq/",m,"_1.fastq",
                  " -2 /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/rawdata/fastq/",m,"_2.fastq",
                  " -S /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/mapping/",m,".sam") %>% 
      as.data.frame()
    df = rbind(df, comm)
  }
  
  comm1 = paste0("samtools sort -@ 70 -o /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/bam/",m,".sorted.bam", 
                 " /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/mapping/",m,".sam") %>% 
    as.data.frame()
  df = rbind(df, comm1)
  
  comm2 = paste0("samtools index /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/bam/",m,".sorted.bam") %>% 
    as.data.frame()
  df = rbind(df, comm2)
  
  comm3 = paste0("htseq-count -f bam -r name -i gene_id -s reverse -t exon -m intersection-nonempty",
                 ' /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/bam/' , m, '.sorted.bam', 
                 ' /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/genome/Nip.gtf',
                 ' > ',  '/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/counts/',  m, '.count.txt') %>% 
    as.data.frame()
  
  if (sample$group[i] == "single") {
    comm4 = paste0("featureCounts -t exon -g gene_id -a /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/genome/Nip.gtf -o ",
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/counts/",m,".counts.featurecounts.txt ",
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/bam/",m,".sorted.bam")  %>% 
      as.data.frame()
    df = rbind(df, comm4)
  }else{
    comm4 = paste0("featureCounts -p -t exon -g gene_id -a /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/data/genome/Nip.gtf -o ",
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/counts/",m,".counts.featurecounts.txt ",
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/bam/",m,".sorted.bam")  %>% 
      as.data.frame()
    df = rbind(df, comm4)
  }
  
  comm5 = paste0("cut -f 1,7 /home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/counts/",m,".counts.featurecounts.txt | grep -v '^#' > ", 
                 "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/counts/",m,".counts.txt") %>% 
    as.data.frame()
  df = rbind(df, comm5)
}


data.table::fwrite(df,
                   "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/stress/mapping.sam2bam.and.featurecounts.2.sh", 
                   col.names = FALSE, row.names = FALSE, quote = FALSE)













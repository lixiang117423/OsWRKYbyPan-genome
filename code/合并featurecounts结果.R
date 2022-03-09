rm(list = ls())

library(GenomicFeatures)
library(data.table)
library(dplyr)

df_all <- data.table::fread("rna-seq/counts/SRR9969498.counts.txt", header = TRUE)
colnames(df_all) <- c("geneid", "SRR9969498")

file_dir <- dir("rna-seq/counts/")

for (i in file_dir) {
  len <- length(stringr::str_split(i, "\\.")[[1]])
  if (len == 3) {
    sample <- stringr::str_split(i, "\\.")[[1]][1]
    df_temp <- data.table::fread(paste0("rna-seq/counts/", i), header = TRUE)
    colnames(df_temp) <- c("geneid", sample)

    if (sample == "SRR9969498") {
      next
    } else {
      df_all <- merge(df_all, df_temp, by = "geneid")
    }
  }
}

data.table::fwrite(df_all,
  file = "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/counts.txt",
  row.names = FALSE, quote = FALSE
)

# counts转换成tpm和FPKM
# txdb <- makeTxDbFromGFF("/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/data/genome/ir64.gff",format="auto")
# xons_gene <- exonsBy(txdb, by = "gene")
# exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
# length=t(as.data.frame(exons_gene_lens))
# write.table(length,'/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/data/genome/ir64.gene.length.txt',col.names=F,row.names=T,quote=F,sep='\t')

gene.len <- fread("/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/data/genome/ir64.gene.length.txt",
  header = FALSE
) %>%
  dplyr::rename(geneid = V1, length = V2)

df.all <- merge(df_all, gene.len, by = "geneid")

kb <- df.all$length / 1000

rpk <- df.all[, 2:(ncol(df.all) - 1)] / kb
rownames(rpk) <- df.all$gene

# TPM
tpm <- t(t(rpk) / colSums(rpk) * 1000000) %>% as.data.frame()
rownames(tpm) <- df.all$gene
data.table::fwrite(tpm, file = "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/TPM.csv", quote = F, row.names = TRUE)

# FPKM
fpkm <- t(t(rpk) / colSums(df_all[, 2:ncol(df_all)]) * 10^6) %>% as.data.frame()
rownames(fpkm) <- df.all$gene
data.table::fwrite(fpkm, file = "/home/lixiang/project/33PanGenomeWRKYFamily/rna-seq/FPKM.csv", quote = F, row.names = TRUE)









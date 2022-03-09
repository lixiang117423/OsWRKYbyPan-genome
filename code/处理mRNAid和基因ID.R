rm(list = ls())


# 提取mRNA2geneID
file_all <- dir("/opt/publicdata/refgenome/rice/33pan_genone/")

for (i in file_all) {
  if (stringr::str_sub(i, 1, 3) != "all") {
    dir_file <- paste0("/opt/publicdata/refgenome/rice/33pan_genone/", i, "/")
    dir_file <- dir(dir_file)
    for (j in dir_file) {
      if (stringr::str_sub(j, nchar(j) - 2, nchar(j)) == "gff") {
        command <- paste0(
          " ../../../code/mRNAid_to_geneid.pl ",
          paste0("/opt/publicdata/refgenome/rice/33pan_genone/", i, "/", j, " "),
          stringr::str_split(j, "\\.")[[1]][1], ".mRNA2geneID.txt"
        )
        system2(command = "perl", args = command)
      }
    }
  }
}

# 查看提取的结果是否正常
files = dir("./")

for (i in files) {
  print(stringr::str_split(i, "\\.")[[1]][1])
  data.table::fread(i, nrows = 5, header = FALSE) %>% 
    print()
  print("----------------------------------------------------------------------")
}


# 确定没有问题后全部合并
final_res = NULL

files = dir("./")

for (i in files) {
  data.table::fread(i, header = TRUE) %>% 
    rbind(final_res) -> final_res
  print("----------------------------------------------------------------------")
}

data.table::fwrite(final_res, 
                   file = "../../results/step1_get_id/mRNA2geneID.txt", 
                   row.names = FALSE,
                   col.names = TRUE,
                   quote = FALSE,
                   sep = "\t"
                   )

################################################################################
# 提取geneID2mRNAid
file_all <- dir("/opt/publicdata/refgenome/rice/33pan_genone/")

for (i in file_all) {
  if (stringr::str_sub(i, 1, 3) != "all") {
    dir_file <- paste0("/opt/publicdata/refgenome/rice/33pan_genone/", i, "/")
    dir_file <- dir(dir_file)
    for (j in dir_file) {
      if (stringr::str_sub(j, nchar(j) - 2, nchar(j)) == "gff") {
        command <- paste0(
          " ../../../code/geneid_to_mRNAid.pl ",
          paste0("/opt/publicdata/refgenome/rice/33pan_genone/", i, "/", j, " "),
          stringr::str_split(j, "\\.")[[1]][1], ".geneID2mRNAid.txt"
        )
        system2(command = "perl", args = command)
      }
    }
  }
}



# 查看提取的结果是否正常
files = dir("./")

for (i in files) {
  print(stringr::str_split(i, "\\.")[[1]][1])
  data.table::fread(i, nrows = 5, header = FALSE) %>% 
    #dplyr::select(1:6) %>% 
    print()
  print("----------------------------------------------------------------------")
}



# 确定没有问题后全部合并
# 合并脚本
head("../../../code/合并mRNA2geneID.py")






















rm(list = ls())

file_dir = dir("/opt/publicdata/refgenome/rice/33pan_genone/")

for (i in file_dir) {
  if (stringr::str_sub(i,1,3) == "all") {
    next
  }else{
    file_dir_2 = dir(paste0("/opt/publicdata/refgenome/rice/33pan_genone/",i))
    for (j in file_dir_2) {
      if (stringr::str_sub(j, nchar(j)-2, nchar(j)) == "gff") {
        print(j)
        
        data.table::fread(paste0("/opt/publicdata/refgenome/rice/33pan_genone/",i,"/",j)) %>% 
          dplyr::mutate(V1 = paste0(V1, "_",i)) %>% 
          data.table::fwrite(file = paste0("data/gff/",i, ".gff"), sep = "\t",
                             row.names = FALSE, col.names = FALSE, quote = FALSE)
      }else{
        next
      }
    }
  }
}


for (i in file_dir) {
  if (stringr::str_sub(i,1,3) == "all") {
    next
  }else{
    file_dir_2 = dir(paste0("/opt/publicdata/refgenome/rice/33pan_genone/",i))
    for (j in file_dir_2) {
      if (stringr::str_sub(j, nchar(j)-1, nchar(j)) == "fa") {
        print(j)
      }else{
        next
      }
    }
  }
}

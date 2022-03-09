import os

f_pep = open("/opt/publicdata/refgenome/rice/33pan_genone/all.pep.fasta","r")

f_new = open("/home/lixiang/project/33PanGenomeWRKYFamily/data/33pangenome/all.pep.new.fasta","w")

f_read = f_pep.readlines()

for line in f_read:
  if " predict " in line:
    line = line.replace(" predict ","")
  else:
    line = line
  
  if "LOC_" in line:
    line = line.split(" ")[0] + "\n"
  else:
    line = line
  
  line = line.replace("P0","T0")
  
  f_new.write(line)
  
f_new.close()


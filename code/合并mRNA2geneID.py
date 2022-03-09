import os

print(os.listdir())
print(os.getcwd())

f = open("../../../results/step1_get_id/geneid2mrnaid.txt","w")

for i in os.listdir():
  if i == ".Rhistory":
    next
  else:
    f_open = open(i, "r")
    f_read = f_open.readlines()
    for line in f_read:
      f.write(line)
      
f.close()


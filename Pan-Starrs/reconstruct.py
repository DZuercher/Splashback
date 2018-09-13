import sys
import frogress

base_name = "/work/dominik.zuercher/DataStore/Pan-Starrs/Chunked_randoms/PS_randoms_21.5"
Outdir = "/work/dominik.zuercher/DataStore/Pan-Starrs"

name_list = [base_name + "/PS_randoms_21.5.dat%03d.dat" % i for i in range(0,980)]
fw = open(Outdir + "/PS_randoms_21.5_reconstructed.dat", "w+")
for name in frogress.bar(name_list):
    with open(name, "r") as f:
	for line in f:
	    fw.write(line)
    f.close()
fw.close()

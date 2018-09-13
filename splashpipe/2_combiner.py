import pandas
import numpy as np
import sys
import glob
from subprocess import call


if __name__ == "__main__":

    dire = sys.argv[1]

    for ii, fname in enumerate(glob.glob(dire+"/???/pairs.dat")):
	print fname
	if ii==0:
	    df = pandas.read_csv(fname, delim_whitespace=1, header=None, names=(["rr", "dd", "avr", "area"]), comment="#")
	    #df_denom = pandas.read_csv(fname+".denom", delim_whitespace=1, header=None, names=(["denom"]), comment="#")
	else:
	    dfnew = pandas.read_csv(fname, delim_whitespace=1, header=None, names=(["rr", "dd", "avr", "area", "denom"]), comment="#")
	    #df_denom_new = pandas.read_csv(fname+".denom", delim_whitespace=1, header=None, names=(["denom"]), comment="#")
	    df.dd = df.dd + dfnew.dd
	    #df_denom.denom = df_denom.denom + df_denom_new.denom
    #df.join(df_denom)
    #df = pandas.concat([df,df_denom],axis=1)
    df.to_csv(dire+"/pairs.dat", index=False, header=False, sep=" ", float_format="%.6e")

    call("cp %s/000/Ncluster_0.dat %s/Ncluster.dat"%(dire, dire), shell=True)

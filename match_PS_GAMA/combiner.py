import pandas
import glob
import numpy as np


if __name__ == "__main__":
    input_dir = "/work/dominik.zuercher/Output/match_PS_GAMA/spec_parts"
    output_dir = "/work/dominik.zuercher/Output/match_PS_GAMA"

    zred=np.zeros(0)
    color_PS=np.zeros(0)

    for ii, fname in enumerate(glob.glob("%s/matched_???" % input_dir)):
        print(fname)
        dfnew = pandas.read_csv(fname, delim_whitespace=1, header=None,names=["zred","color_PS"])
        zred=np.append(zred,dfnew.zred.values)
        color_PS = np.append(color_PS,dfnew.color_PS.values)
    np.savetxt("%s/matched_spec.dat" % output_dir, np.hstack((zred.reshape(zred.size,1),color_PS.reshape(color_PS.size,1))) )


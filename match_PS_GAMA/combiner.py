import pandas
import glob
import numpy as np


if __name__ == "__main__":
    input_dir = "/work/dominik.zuercher/Output/match_PS_GAMA/spec_parts"
    output_dir = "/work/dominik.zuercher/Output/match_PS_GAMA"

    zred=np.zeros(0)
    red_PS=np.zeros(0)
    green_PS=np.zeros(0)
    iband_PS=np.zeros(0)

    for ii, fname in enumerate(glob.glob("%s/matched_???" % input_dir)):
        print(fname)
        dfnew = pandas.read_csv(fname, delim_whitespace=1, header=None,names=["zred","red_PS", "green_PS", "iband_PS"])
        zred=np.append(zred,dfnew.zred.values)
        red_PS=np.append(red_PS,dfnew.red_PS.values)
        green_PS=np.append(green_PS,dfnew.green_PS.values)
        iband_PS=np.append(iband_PS,dfnew.iband_PS.values)
    np.savetxt("%s/matched_spec.dat" % output_dir, np.hstack((zred.reshape(zred.size,1), red_PS.reshape(red_PS.size,1), green_PS.reshape(green_PS.size,1), iband_PS.reshape(iband_PS.size,1) )) )


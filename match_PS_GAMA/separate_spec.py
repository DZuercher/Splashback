import pandas as pd
import numpy as np


if __name__ == "__main__":

    dir_ = "/work/dominik.zuercher/Output/match_PS"

    spec_cat = "%s/matched_spec.dat" % dir_
    data = pd.read_csv(spec_cat, delim_whitespace=1,header=None,names=["zred","red","green","red_PS","green_PS","iband_PS"])
    cat = data.values
    id_ = np.zeros(cat.shape[0])
    id_ = (cat[:,2] - cat[:,1] <= 0.8 - 0.03*(cat[:,1] + 20.0)) #True if identified as blue and False if red
    cat = np.hstack((cat,id_.reshape(id_.size,1))) 

    id_ = np.zeros(cat.shape[0])
    id_ = (cat[:,4] - cat[:,3] <= 0.8 - 0.03*(cat[:,3] + 20.0)) #True if identified as blue and False if red
    cat = np.hstack((cat,id_.reshape(id_.size,1))) 
    np.savetxt("%s/matched_spec_new.dat" % dir_,cat)


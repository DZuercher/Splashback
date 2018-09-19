import numpy as np
import pandas as pd

def combine(type_, modded=False):

    if modded==False:
        value_array=pd.read_csv("%s/%s/Data_full.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)
    else:
        value_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    value_array = np.asarray(value_array.values)


    dev2_array = pd.read_csv("%s/%s/dev2_full.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    dev2_array = np.asarray(dev2_array.values)

    np.hstack((value_array, dev2_array))

   
    np.savetxt("%s/%s/Data_full_modded_dev2.txt" % (output_dir, type_), np.hstack((value_array, dev2_array)) )


if __name__ == "__main__":

    output_dir = "/work/dominik.zuercher/Output/Mest"
    types=["Planck_PS_21.5_red_spline", "Planck_PS_21.5_blue_spline"]
    adds = ["_best"]
    modded = True
    for type_ in types:
	combine(type_,modded)


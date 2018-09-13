import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import frogress


def correct_chain(type_,add,modded=False):
    #Read calculated values
    print("Reading MCMC_data...")
    
    if modded==False:
        value_array=pd.read_csv("%s/%s/Data_full.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)
    else:
        value_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    value_array=np.asarray(value_array.values)
    #Read chains
    print("Reading chains...")
    if modded==False:
        samples=np.loadtxt("%s/%s/chainstate_full.txt" % (output_dir, type_))

    else:
        data = pd.read_csv("%s/%s/chainstate_full_modded.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    samples=samples.astype(float)
    #Find broken chains (do visual inspection of chain file and come up with criterion)
    print("Removing broken chains...")
    if broken_walker!=None:
	broken_walkers=[broken_walker+i*nwalkers for i in range(steps)]
	#Remove broken chains
	value_array=np.delete(value_array,broken_walkers,axis=0)
	samples=np.delete(samples,broken_walkers,axis=0)

    #Find nan values from rhodev_array
    print("Search and remove NAN entries...")
    rhodev_array=value_array[:,rsteps*3:rsteps*4]
    fix=np.isnan(rhodev_array)
    row=np.invert(np.any(fix,axis=1))
    print("Found %d NAN derivatives" % (value_array.shape[0] - np.sum(row)) )

    #Remove nan value entries from rhodev_array
    value_array=value_array[row,:]

    print("Saving modded MCMC_data...")
    np.savetxt("%s/%s/Data_full_modded.txt" % (output_dir, type_), value_array)

    samples=samples.astype(float)
    print("Saving modded chains...")
    np.savetxt("%s/%s/chainstate_full_modded.txt" % (output_dir, type_), samples)


if __name__=="__main__":

    output_dir = "/work/dominik.zuercher/Output/Mest"

    steps = 11000
    nwalkers = 28
    rsteps = 25
    broken_walker = None   


    types=["Planck_PS_21.5_red_spline", "Planck_PS_21.5_blue_spline"]
    adds = ["_best"]
    modded = False
    for add in adds:
        for type_ in types:
            correct_chain(type_,add,modded)


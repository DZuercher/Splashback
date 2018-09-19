import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/Splashback/toolbox')
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
    if mc == False:
        type_ += "_no_mc"

    #Read calculated values
    print("Reading MCMC_data...")
    
    if modded==False:
        value_array=pd.read_csv("%s/%s/Data_full.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)
    else:
        value_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    value_array=np.asarray(value_array.values)


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


    samples=samples.astype(float)
    print("Saving modded chains...")
    np.savetxt("%s/%s/chainstate_full_modded.txt" % (output_dir, type_), samples)

def correct_dev2(type_,modded):
    if mc == False:
        type_ += "_no_mc"

    if modded==False:
        dev2_array=pd.read_csv("%s/%s/dev2_full.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)
    else:
        dev2_array=pd.read_csv("%s/%s/dev2_full_modded.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    dev2_array=np.asarray(dev2_array.values)
    #Find nan values from rhodev2_array
    print("Search and remove NAN entries...")
    fix=np.isnan(dev2_array)
    row=np.invert(np.any(fix,axis=1))
    print("Found %d NAN derivatives" % (dev2_array.shape[0] - np.sum(row)) )

    #Remove nan value entries from rhodev2_array
    dev2_array=dev2_array[row,:]
    print("Saving modded MCMC_data...")
    np.savetxt("%s/%s/dev2_full_modded.txt" % (output_dir, type_), dev2_array)

if __name__=="__main__":

    output_dir = "/work/dominik.zuercher/Output/Mest"

    steps = 100000
    nwalkers = 28
    rsteps = 25
    broken_walker = None   

    mc = False

    types=["Planck_PS_21.5_red_spline", "Planck_PS_21.5_blue_spline", "Planck_PS_21", "Planck_PS_21.5", "Planck_PS_22"]
    adds = ["_best"]
    modded = False
    for add in adds:
        for type_ in types:
            #correct_chain(type_,add,modded)
	    correct_dev2(type_,modded)

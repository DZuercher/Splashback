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



def cal_procedure(type_,add,modded=False):

    #Reading data
    if modded==False:
        data_array=pd.read_csv("%s/%s/Data_full.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)
    else:
        data_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (output_dir, type_), sep=' ',header=None,error_bad_lines=False)

    print("Read")
    data_array=np.asarray(data_array.values)
    sigma_array=data_array[:,0:rsteps]
    rho_array=data_array[:,rsteps:rsteps*2]
    dev_array=data_array[:,rsteps*2:rsteps*3]
    rhodev_array=data_array[:,rsteps*3:rsteps*4]
    rsp2d_array=data_array[:,rsteps*4]
    rsp2d_array=rsp2d_array.reshape((rsp2d_array.size,1))
    rsp3d_array=data_array[:,rsteps*4+1]
    rsp3d_array=rsp3d_array.reshape((rsp3d_array.size,1))
    chisquare_array=data_array[:,rsteps*4+2]
    chisquare_array=chisquare_array.reshape((chisquare_array.size,1))
    maxrhodev_array=data_array[:,rsteps*4+3]
    maxrhodev_array=maxrhodev_array.reshape((maxrhodev_array.size,1))
    maxinnerrhodev_array=data_array[:,rsteps*4+4]
    maxinnerrhodev_array=maxinnerrhodev_array.reshape((maxinnerrhodev_array.size,1))


    #Calculating percentile values
    sigmas = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(sigma_array, [16, 50, 84], axis=0))))
    print("sigmas done")
    rhos = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(rho_array, [16, 50, 84], axis=0))))
    print("rhos done")
    devs = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(dev_array, [16, 50, 84], axis=0))))
    print("devs done")
    rhodevs = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(rhodev_array, [16, 50, 84], axis=0))))
    print("rhodevs done")
    rsp2ds = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(rsp2d_array, [16, 50, 84], axis=0))))
    rsp2ds=rsp2ds.reshape(rsp2ds.size)
    print("rsp2ds done")
    rsp3ds = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(rsp3d_array, [16, 50, 84], axis=0))))
    rsp3ds=rsp3ds.reshape(rsp3ds.size)
    print("rsp3ds done")
    chisquares = np.min(chisquare_array)
    print(chisquares)
    chisquares= np.asarray([chisquares,0,0])
    chisquares = chisquares.reshape(chisquares.size)
    print("chisquares done")
    maxrhodevs = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(maxrhodev_array, [16, 50, 84], axis=0))))
    maxrhodevs=maxrhodevs.reshape(maxrhodevs.size)
    print("maxrhodevs done")
    maxinnerrhodevs = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(maxinnerrhodev_array, [16, 50, 84], axis=0))))
    maxinnerrhodevs=maxinnerrhodevs.reshape(maxinnerrhodevs.size)
    print("maxinnerrhodevs done")

    #Saving
    outarr=np.vstack((sigmas,rhos,devs,rhodevs))
    np.savetxt("%s/%s/plotsave.dat" % (output_dir, type_), outarr)
    valarr=np.vstack((rsp2ds,rsp3ds,chisquares,maxrhodevs,maxinnerrhodevs))
    np.savetxt("%s/%s/values.dat" % (output_dir, type_), valarr)


    #Plotting Histograms
    f, axarr = plt.subplots(2,2)
    n,bins,patches = axarr[0,0].hist(rsp2d_array.reshape(rsp2d_array.size),100,alpha=0.5,range=[0,10],label="RSP 2D")
    axarr[0,0].hist(rsp3d_array.reshape(rsp3d_array.size),bins,alpha=0.5,label="RSP 3D")
    axarr[0,0].legend()
    axarr[0,0].set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    axarr[0,0].set_xlim([0,6])

    n,bins,patches = axarr[1,0].hist(maxrhodev_array.reshape(maxrhodev_array.size),50,range=[-5,0],alpha=0.5,label="Max Rho dev")
    nans=np.invert(np.isnan(maxinnerrhodev_array))
    maxinnerrhodev_array=maxinnerrhodev_array[nans.reshape((nans.size))]
    axarr[1,0].legend()

    n,bins,patches=axarr[0,1].hist(chisquare_array.reshape(chisquare_array.size),bins=500,range=[0,500],alpha=0.5,label="CHI2")
    axarr[0,1].legend()

    plt.tight_layout()
    plt.savefig("%s/%s/vis_values.pdf" % (output_dir, type_))



#-----Main-----

if __name__ == "__main__":

    rsteps = 25
    mc = True
    modded = True

    output_dir = "/work/dominik.zuercher/Output/Mest"

    types = ["Planck_PS_21.5_blue_spline","Planck_PS_21.5_red_spline"]
    adds = ["_best"]
    for add in adds:
        for type_ in types:
            print("Doing "+str(type_))
            cal_procedure(type_,add,modded)

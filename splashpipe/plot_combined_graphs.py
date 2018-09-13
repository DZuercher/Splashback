import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def compare(ests,types,name,pres):
    colors=np.linspace(0,1,int(len(ests)*len(types)))
    if pres==True:
        #shifts=np.linspace(1,10+(int(len(ests)*len(types))-1),int(len(ests)*len(types)))
        shifts=[9]*(int(len(ests)*len(types)))
        #shifts[0]+=10
    else:
        shifts=np.asarray([1]*int(len(ests)*len(types)))
    fig=plt.figure(1)
    ax = plt.subplot(111)
    i=0
    for type_ in types:
        for est in ests:
            if type_=='surhud':
                dat=np.genfromtxt("/home/dominik.zuercher/Documents/RSP_Pro/Mest/redmap.dat")
            elif "SDSS" in type_:
                     dat=np.genfromtxt("/work/dominik.zuercher/DataStore/corr-pairs/Planck_SDSS/Planck_SDSS_plot("+str(est)+").dat")
            else:
                 try:
                     dat=np.genfromtxt("/work/dominik.zuercher/DataStore/corr-pairs/"+str(type_)+"/"+str(type_)+"_plot("+str(est)+").dat")
                 except:
                     continue
            print("-----------------------------------------")
            print(type_,est)
            print(dat)
            print("-----------------------------------------")
            if type_=='surhud':
                ax.errorbar(dat[:,0], np.multiply(dat[:,1],0.1), np.multiply(dat[:,2],0.1), fmt=".",c=cm.nipy_spectral(colors[i]),capsize=2, alpha=1,label ="RedMaPPer") 
            else:
                ax.errorbar(dat[:,0], np.multiply(dat[:,1],shifts[i]), np.multiply(dat[:,2],shifts[i]), fmt=".",capsize=2,c=cm.nipy_spectral(colors[i]), alpha=0.8,label =str(type_)+" ("+str(est)+")")
            i+=1
    ax.set_title("Comparison")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    ax.set_ylabel(r"$\xi^{\rm 2d}$ ($h^{-1}$Mpc)")
    plt.legend()
    plt.savefig("comparisons/"+str(name)+".pdf")
    plt.close()

def sigma_compare(ests,types,name):
    colors=np.linspace(0,1,int(len(ests)*len(types)))
    fig=plt.figure(2)
    ax = plt.subplot(111)
    i=0
    for type_ in types:
        for est in ests:
            if type_=='surhud':
                dat=np.genfromtxt("/home/dominik.zuercher/Documents/RSP_Pro/Mest/redmap.dat")
            elif "SDSS" in type_:
                dat=np.genfromtxt("/work/dominik.zuercher/DataStore/corr-pairs/Planck_SDSS/Planck_SDSS_plot("+str(est)+").dat")
            else:
                try:
                    dat=np.genfromtxt("/work/dominik.zuercher/DataStore/corr-pairs/"+str(type_)+"/"+str(type_)+"_sigplot("+str(est)+").dat")
                except:
                    continue
            print("-----------------------------------------")
            print(name)
            print(dat)
            print("-----------------------------------------")
            if type_=='surhud':
                ax.errorbar(dat[:,0], dat[:,1], dat[:,2], fmt=".",c=cm.nipy_spectral(colors[i]),capsize=2, alpha=1,label ="RedMaPPer")
            else:
                ax.errorbar(dat[:,0], dat[:,1], dat[:,2], fmt=".",c=cm.nipy_spectral(colors[i]),capsize=2, alpha=0.8,label =str(type_)+" ("+str(est)+")")
            i+=1
    ax.set_title("Sigma Comparison")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R$ ($h^{-1}$Mpc)")
    ax.set_ylabel(r"$\Sigma_g$ ($h^{2}$Mpc^{-2})")
    plt.legend()
    plt.savefig("sigmacomparisons/"+str(name)+".pdf")
    plt.close()



#types=["Planck_PS",'Planck_PS_area','Planck_PS_full','Planck_PS_full_area','XR_PS','XR_PS_area','SDSS']
#types=['surhud','Planck_PS_21','Planck_PS_21.5','Planck_PS_22']
#ests=['DR'] #,'RD','RR']

pres=True
"""
#Depth comparioson plots
for type_ in ["Planck_PS",'Planck_PS_area','Planck_SDSS']:
   for est in ["DR","RD","RR","R","D"]:
       nums=["_21","_21.5","_22"]
       types_=[type_+num for num in nums]
       name=type_+"_"+est
       print("Doing "+str(name))
       est=[est]
       compare(est,types_,name,pres)
       #sigma_compare(est,types_,name)

#est comparion plots
for type_ in ["Planck_PS",'Planck_PS_area','Planck_SDSS']:
   for num in ["_21","_21.5","_22"]:
       ests=["DR","RD","RR","R","D"]
       types_=type_+num
       types_=[types_]*3
       name=type_+num
       print("Doing "+str(name))
       compare(ests,types_,name,pres)
       #sigma_compare(ests,types_,name)

#Type comparison plots
for est in ["DR","RD","RR","R","D"]:
    est=[est]
    types=["Planck_PS",'Planck_PS_area','Planck_SDSS']
    for num in ["_21","_21.5","_22"]:
        types_=[type_+num for type_ in types]
        name=est[0]+num
        print("Doing "+str(name))
        compare(est,types_,name,pres)
        #sigma_compare(est,types_,name)
"""
#------------------------------
#Report plots
#------------------------------

est=["DR","RR"]
types=["Planck_PS_21_rebin","Planck_PS_21.5_rebin","Planck_PS_22_rebin","surhud"]
compare(est,types,name="presentation",pres=True)


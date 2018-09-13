import numpy as np
import pandas as pd
import scipy.linalg
def calc_area(type_):
   if type_=='Planck_PS_21':
       randinit=1304704002
       randafter=617205099
   elif type_=='Planck_PS_22':
       randinit=1155286477.0*4.0/3.0
       randafter=585841750.0
   elif type_=='Planck_PS_21.5':
       randinit=1155286477.0*4.0/3.0
       randafter=768660433
   area=randafter/(randinit/(41253))
   return area


def calc_SNR(type_):
    if type_=='surhud':
        data=np.genfromtxt("/home/dominik.zuercher/Documents/RSP_Pro/Mest/redmap.dat")
        cov_data = np.genfromtxt("/home/dominik.zuercher/Documents/RSP_Pro/Mest/redmap.cov")[:,0]
        cov=cov_data.reshape((data[:,0].size,data[:,0].size))
    else:
	data = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro_reload/splashpipe/reloaded/"+str(type_)+"/xi_2d.dat", header=None, sep = ' ')
	cov_data = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro_reload/splashpipe/reloaded/"+str(type_)+"/xi_2d_cov.dat", header=None, sep = ' ')
        cov = np.asarray(cov_data.values)
        data = np.asarray(data.values)

    rr = data[:,0]
    yy = data[:,1]
    yy_err = data[:,2]

    invcov=scipy.linalg.inv(cov)

    SNR=np.sqrt(np.dot(yy,np.dot(invcov,yy)))
    #SNR=np.sqrt(np.dot(yy.T, np.dot(invcov, yy)))
    return SNR



#Survey areas
types=['Planck_PS_21','Planck_PS_22','Planck_PS_21.5']
for type_ in types:
    print("Effective survey area of "+str(type_)+" is "+str(calc_area(type_)))

#SNRs
#types=["Planck_PS_21",'Planck_PS_area_21','Planck_PS_full_21','Planck_PS_full_area_21','XR_PS_21','XR_PS_area_21',"Planck_PS_21.5",'Planck_PS_area_21.5','Planck_PS_full_21.5','Planck_PS_full_area_21.5','XR_PS_21.5','XR_PS_area_21.5',"Planck_PS_22",'Planck_PS_area_22','Planck_PS_full_22','Planck_PS_full_area_22','XR_PS_22','XR_PS_area_22','Planck_SDSS','surhud']
types=["Planck_PS_21_reloaded","Planck_PS_21.5_reloaded","Planck_PS_22_reloaded",'surhud']
for type_ in types:
    SNR=calc_SNR(type_)
    if SNR==1:
	pass
    else:
	print("The SNR for "+str(type_)+" is "+str(SNR))







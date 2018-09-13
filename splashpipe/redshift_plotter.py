import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
dirin="/work/dominik.zuercher/DataStore/corr-pairs/Planck_PS_21_rebin/Normal"
full=np.zeros(0)
for i in range(0,30):
	normal=np.loadtxt(dirin+"/redshift_distribution_"+str(i)+".txt")
        print(normal)
        full=np.append(full,normal)
        """
	for rank in range(0,168):
		data=np.load(dirin+"/save_random_"+str(rank)+"_"+str(i)+".npy")
		fulldata=np.add(fulldata,data[0,:])
	fulldatasum=fulldata[-1]
	fulldata=np.divide(fulldata,fulldatasum)
        """
plt.hist(full)
plt.savefig(dirin+"/redshift_dist.pdf")

# Splits the cluster catalog into jackknife bins based on their locations on the sky
# Each bin contains a similiar number of clusters

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def fullsplit(vals, regions):
    splits = np.asarray(np.array_split(vals, regions))
    middles = np.zeros(0)
    for it in range(splits.shape[0] - 1):
        middle = np.median(np.append(splits[it], splits[it+1]))
	middles = np.append(middles, middle)
    return np.append(np.append(-40.0, middles), 360.0), splits


def get_north_lines(npoints):
    nras = npoints[:,0]
    ndecs = npoints[:,1]

    sortid = np.argsort(nras)
    nras = np.take(nras, sortid)
    ndecs = np.take(ndecs, sortid)

    nralines, nrasplits = fullsplit(nras, 7)
    totlength = 0
    argdecsplits = np.zeros(7, dtype = object)
    for i in range(nrasplits.shape[0]):
	length = nrasplits[i].size
        new = length + totlength
        argdecsplits[i] = np.sort(ndecs[totlength:new])
        totlength = new

    north = np.zeros((8 + 1, 1))
    for i, part in enumerate(argdecsplits):
        declines, decsplits = fullsplit(part, 8)
        declines = declines.reshape((declines.size, 1))
        north = np.hstack((north, declines))
    north = north[:,1:]
    return nralines, north


def get_south_lines(spoints):
    sras = spoints[:,0]
    #shift to patch southern part together
    sras = np.mod(np.add(sras, 100.0), 360.0)
    sdecs = spoints[:,1]

    sortid = np.argsort(sras)
    sras = np.take(sras, sortid)
    sdecs = np.take(sdecs, sortid)

    sralines, srasplits = fullsplit(sras, 6)

    totlength = 0
    argdecsplits = np.zeros(6, dtype = object)
    for i in range(srasplits.shape[0]):
        length = srasplits[i].size
        new = length + totlength
        argdecsplits[i] = np.sort(sdecs[totlength:new])
        totlength = new

    south = np.zeros((7 + 1, 1))
    for i, part in enumerate(argdecsplits):
        declines, decsplits = fullsplit(part, 7)
        declines = declines.reshape((declines.size, 1))
        south = np.hstack((south, declines))
    south = south[:,1:]
    return sralines, south


def northsplit(ncl_ras, ncl_decs):
    jackregs = np.zeros_like((ncl_ras))
    njackreg = 56 #regions in northern part
    rasplit = 7
    decsplit = 8
    for it in range(njackreg):
        decit = int(it/rasplit)
        rait = it%rasplit
        cl_ras_reg = (ncl_ras >= nralines[rait]) & (ncl_ras < nralines[rait + 1])
        cl_decs_reg = (ncl_decs >= north[decit,rait]) & (ncl_decs < north[decit + 1, rait])
        idx = (cl_ras_reg==1) & (cl_decs_reg==1)
        jackregs[idx] = int(it)
    return jackregs


def southsplit(scl_ras, scl_decs):
    scl_ras = np.mod(np.add(scl_ras, 100.0), 360.0)
    jackregs = np.zeros_like((scl_ras))
    sjackreg = 42 #regions in southern part
    rasplit = 6
    decsplit = 7
    for it in range(sjackreg):
        decit = int(it/rasplit)
        rait = it%rasplit
        cl_ras_reg = (scl_ras >= sralines[rait]) & (scl_ras < sralines[rait + 1])
        cl_decs_reg = (scl_decs >= south[decit,rait]) & (scl_decs < south[decit + 1, rait])
        idx = (cl_ras_reg==1) & (cl_decs_reg==1)
        jackregs[idx] = int(it)
    return np.add(jackregs, 56)


if __name__=="__main__":

print("-----------------Splitting PS randoms---------------")
print("Splitting North cap")
#North cap is 7X8 division
npoints = pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/PS_randoms_north.dat",header=None,sep=" ")
npoints=np.asarray(npoints.values)

nralines, north=get_north_lines(npoints)

np.savetxt("north_ras.txt",nralines)
np.savetxt("north_decs.txt",north)


print("Splitting South cap")
#South cap is 6X7 division
spoints=pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/PS_randoms_south.dat",header=None,sep=" ")
spoints=np.asarray(spoints.values)

sralines, south=get_south_lines(spoints)

np.savetxt("south_ras.txt",sralines)
np.savetxt("south_decs.txt",south)

print("----------Northern cuts---------")
print(nralines)
print(north)
print("----------Southern cuts---------")
print(sralines)
print(south)













print("------------------Splitting Planck clusters---------------------")
ncl=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_north.csv",sep=" ")
scl=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_south.csv",sep=" ")


print("Splitting North cap")
pncl_ras=np.asarray(ncl["RA"].values)
pncl_decs=np.asarray(ncl["DEC"].values)

jackregs=northsplit(pncl_ras,pncl_decs)
ncl["jackreg"]=jackregs


print("Splitting South cap")
pscl_ras=np.asarray(scl["RA"].values)
pscl_decs=np.asarray(scl["DEC"].values)

jackregs=southsplit(pscl_ras,pscl_decs)
scl["jackreg"]=jackregs


#Combine north and south and save
output=ncl.append(scl)
output.to_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_area.csv",index=False)












"""
print("-----------------Splitting Full Planck clusters-----------------")

ncl=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Full_Planck_north.csv",sep=" ")
scl=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Full_Planck_south.csv",sep=" ")


print("Splitting North cap")
fpncl_ras=np.asarray(ncl["RA"].values)
fpncl_decs=np.asarray(ncl["DEC"].values)

jackregs=northsplit(fpncl_ras,fpncl_decs)
ncl["jackreg"]=jackregs


print("Splitting South cap")
fpscl_ras=np.asarray(scl["RA"].values)
fpscl_decs=np.asarray(scl["DEC"].values)

jackregs=southsplit(fpscl_ras,fpscl_decs)
scl["jackreg"]=jackregs


#Combine north and south and save
output=ncl.append(scl)
output.to_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Full_Planck_area.csv",index=False)
"""

print("-----------------Splitting Random Planck clusters-----------------")

ncl=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_randoms_north.csv",sep=" ")
scl=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_randoms_south.csv",sep=" ")


print("Splitting North cap")
prncl_ras=np.asarray(ncl["RA"].values)
prncl_decs=np.asarray(ncl["DEC"].values)

jackregs=northsplit(prncl_ras,prncl_decs)
ncl["jackreg"]=jackregs


print("Splitting South cap")
prscl_ras=np.asarray(scl["RA"].values)
prscl_decs=np.asarray(scl["DEC"].values)

jackregs=southsplit(prscl_ras,prscl_decs)
scl["jackreg"]=jackregs


#Combine north and south and save
output=ncl.append(scl)
output.to_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_randoms_area.csv",index=False)


print("-----------------Splitting XR clusters-----------------")

ncl=pd.read_csv("/work/dominik.zuercher/DataStore/XR/north_cap.dat",sep=" ")
scl=pd.read_csv("/work/dominik.zuercher/DataStore/XR/south_cap.dat",sep=" ")

print("Splitting North cap")
xncl_ras=np.asarray(ncl["RA"].values)
xncl_decs=np.asarray(ncl["DEC"].values)

jackregs=northsplit(xncl_ras,xncl_decs)
ncl["jackreg"]=jackregs


print("Splitting South cap")
xscl_ras=np.asarray(scl["RA"].values)
xscl_decs=np.asarray(scl["DEC"].values)

jackregs=southsplit(xscl_ras,xscl_decs)
scl["jackreg"]=jackregs


#Combine north and south and save
output=ncl.append(scl)
output.to_csv("/work/dominik.zuercher/DataStore/XR/mcxc_area.csv",index=False)



print("---------------------Plotting------------------")
"""
psnpoints=pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/PS_randoms_north.dat",sep=" ")
psnpoints=psnpoints.values
psnpoints_ra=psnpoints[:,0]
psnpoints_dec=psnpoints[:,1]
psspoints=pd.read_csv("/work/dominik.zuercher/DataStore/Pan-Starrs/Randoms/PS_randoms_south.dat",sep=" ")
psspoints=psspoints.values
psspoints_ra=psspoints[:,0]
psspoints_dec=psspoints[:,1]
"""
f=plt.figure(1)
#Plot random PS points
#plt.scatter(psnpoints_ra,psnpoints_dec,c='b',s=0.1)
#plt.scatter(psspoints_ra,psspoints_dec,c='r',s=0.1)

plt.scatter(pncl_ras,pncl_decs,color="b",s=0.1)
plt.scatter(pscl_ras,pscl_decs,color="r",s=0.1)

#plt.scatter(fpncl_ras,fpncl_decs,color="b",s=0.1)
#plt.scatter(fpscl_ras,fpscl_decs,color="r",s=0.1)

#plt.scatter(prncl_ras,prncl_decs,color="b",s=0.1)
#plt.scatter(prscl_ras,prscl_decs,color="r",s=0.1)

#plt.scatter(xncl_ras,xncl_decs,color="b",s=0.1)
#plt.scatter(xscl_ras,xscl_decs,color="r",s=0.1)

#Plot area borders
#Northern lines
for i in range(nralines.size-1):
	plt.axvline(x=nralines[i],c='k',lw=0.5)
	for j in range(north.shape[0]):
		plt.axhline(y=north[j,i],xmin=nralines[i]/360.0,xmax=nralines[i+1]/360.0,c='k',lw=0.5)
plt.axvline(x=nralines[-1],c='k',lw=0.5)

plotsralines=np.zeros(0)
id0=np.sort(sralines[sralines>=100.0])[0]
id0= np.where(np.abs(sralines-id0)<0.0001)[0]
for i,line in enumerate(sralines):
	if line>=100.0:
		plotsralines=np.append(plotsralines,line-100.0)
	else:
		plotsralines=np.append(plotsralines,line-100.0+360.0)

id_=np.argsort(plotsralines[1:-1])
plotsralines=np.take(plotsralines[1:-1],id_)
plotsralines=np.append(np.append(-40.0,plotsralines),360.0)
southfoo=np.hstack((south[:,:int(id0)],south[:,int(id0)+1:]))
southfoo2=np.take(southfoo,id_,axis=1)
south=np.hstack((south[:,int(id0)].reshape((south.shape[0],1)),southfoo2))

for i in range(plotsralines.size-1):
	plt.axvline(x=plotsralines[i],c='y',lw=0.5)
	for j in range(south.shape[0]):
		plt.axhline(y=south[j,i],xmin=plotsralines[i]/360.0,xmax=plotsralines[i+1]/360.0,c='y',lw=0.5)
plt.axvline(x=plotsralines[-1],c='y',lw=0.5)

plt.ylim((-35.0,90.0))
plt.xlim((0.0,360.0))
plt.savefig("areasplit.pdf")


fig=plt.figure(2)
dat1=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_area.csv",sep=',')
#dat2=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Full_Planck_area.csv",sep=',')
dat3=pd.read_csv("/work/dominik.zuercher/DataStore/Planck-SZ/98/Planck_randoms_area.csv",sep=',')
#dat4=pd.read_csv("/work/dominik.zuercher/DataStore/XR/mcxc_area.csv",sep=',')

colors=np.linspace(0,1,98)
for i in range(98):
	ras=dat1["RA"][dat1["jackreg"]==i].values
	decs=dat1["DEC"][dat1["jackreg"]==i].values
	plt.scatter(ras,decs,c=cm.nipy_spectral(colors[i]),s=0.1)
        """
        ras=dat2["RA"][dat2["jackreg"]==i].values
        decs=dat2["DEC"][dat2["jackreg"]==i].values
        plt.scatter(ras,decs,c=cm.nipy_spectral(colors[i]),s=0.1)

        ras=dat3["RA"][dat3["jackreg"]==i].values
        decs=dat3["DEC"][dat3["jackreg"]==i].values
        plt.scatter(ras,decs,c=cm.nipy_spectral(colors[i]),s=0.1)

        ras=dat4["RA"][dat4["jackreg"]==i].values
        decs=dat4["DEC"][dat4["jackreg"]==i].values
        plt.scatter(ras,decs,c=cm.nipy_spectral(colors[i]),s=0.1)
        """
#Plot area borders
#Northern lines
for i in range(nralines.size-1):
        plt.axvline(x=nralines[i],c='k',lw=0.5)
        for j in range(north.shape[0]):
                plt.axhline(y=north[j,i],xmin=nralines[i]/360.0,xmax=nralines[i+1]/360.0,c='k',lw=0.5)
plt.axvline(x=nralines[-1],c='k',lw=0.5)


for i in range(plotsralines.size-1):
        plt.axvline(x=plotsralines[i],c='y',lw=0.5)
        for j in range(south.shape[0]):
                plt.axhline(y=south[j,i],xmin=plotsralines[i]/360.0,xmax=plotsralines[i+1]/360.0,c='y',lw=0.5)
plt.axvline(x=plotsralines[-1],c='y',lw=0.5)

plt.ylim((-35.0,90.0))
plt.xlim((0.0,360.0))
plt.savefig("jackknife_samples.pdf")

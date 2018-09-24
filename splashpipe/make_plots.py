import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import palettable 
cm=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors


def assign_color(i,idx):
    return cm[i]

def assign_color_xr(i):
    return cm[i*2]

def curve_graphs(types):
    f1, axarr1 = plt.subplots(ncols=2,nrows=3,figsize=(6,6))
    f2, axarr2 = plt.subplots(ncols=3,nrows=3,figsize=(9,6))
    i=0
    shifts=np.asarray([1,1,1])
    for type_ in types:
	print("Doing %s",type_)
	dr_dat = np.genfromtxt("%s/%s/xi_2d.dat" % (splash_directory, type_))
        
        if mc == False:
            type_ += "_no_mc"

	idx = (dr_dat[:,0] > 0.1) & (dr_dat[:,0] < 10.0)
	rr=dr_dat[idx,0]
	rrange=np.linspace(rr[0],rr[-1],rsteps)

	dr_plotdat=pd.read_csv("%s/%s/plotsave.dat" % (est_directory, type_),sep=' ',header=None)
	dr_plotdat=np.asarray(dr_plotdat.values)
	dr_sigmas=dr_plotdat[0:rsteps,:]
	dr_rhos=dr_plotdat[rsteps:rsteps*2,:]
	dr_devs=dr_plotdat[rsteps*2:rsteps*3,:]
	dr_rhodevs=dr_plotdat[rsteps*3:rsteps*4,:]
	dr_rhodevs2=dr_plotdat[rsteps*4:rsteps*5,:]

	dr_splashdat=pd.read_csv("%s/%s/values.dat" % (est_directory, type_), sep=' ',header=None)
	dr_splashdat=np.asarray(dr_splashdat.values)
	dr_rsp2ds=dr_splashdat[0,:]
	dr_rsp3ds=dr_splashdat[1,:]
        
        print("Data read")
        print("Plotting...")

	#Plot 2D data
	axarr1[i,0].errorbar(dr_dat[:,0], dr_dat[:,1], dr_dat[:,2], fmt=".",capsize=2,c='k')

	#Plot 2D MCMC fits
	axarr1[i,0].plot(rrange,dr_sigmas[:,0],linestyle='-',c=assign_color(i,0))
        axarr1[i,0].fill_between(rrange,dr_sigmas[:,0]-dr_sigmas[:,2],dr_sigmas[:,0]+dr_sigmas[:,1],facecolor=assign_color(i,0),alpha=0.3)

	#Make vertical line at locations of R200m and R2D
	axarr1[i,0].axvline(x=r200m,c="k",linestyle=":")
        axarr1[i,0].axvspan(dr_rsp2ds[0]-dr_rsp2ds[2],dr_rsp2ds[0]+dr_rsp2ds[1],color=assign_color(i,0),alpha=0.3)

	axarr1[i,1].axvline(x=r200m,c="k",linestyle=":")
        axarr1[i,1].axvspan(dr_rsp2ds[0]-dr_rsp2ds[2],dr_rsp2ds[0]+dr_rsp2ds[1],color=assign_color(i,0),alpha=0.3)

	#Plot 2D derivative MCMC fits
	axarr1[i,1].plot(rrange,dr_devs[:,0],linestyle='-',c=assign_color(i,0))
        axarr1[i,1].fill_between(rrange,dr_devs[:,0]-dr_devs[:,2],dr_devs[:,0]+dr_devs[:,1],facecolor=assign_color(i,0),alpha=0.3)

	#Plot 3D MCMC fits
        axarr2[i,0].plot(rrange,dr_rhos[:,0],linestyle='-',c=assign_color(i,0))
        axarr2[i,0].fill_between(rrange,dr_rhos[:,0]-dr_rhos[:,2],dr_rhos[:,0]+dr_rhos[:,1],facecolor=assign_color(i,0),alpha=0.3)

	#Make vertical line at locations of R200m and R3D
	axarr2[i,0].axvline(x=r200m,c="k",linestyle=":")
        axarr2[i,0].axvspan(dr_rsp3ds[0]-dr_rsp3ds[2],dr_rsp3ds[0]+dr_rsp3ds[1],color=assign_color(i,0),alpha=0.3)

	axarr2[i,1].axvline(x=r200m,c="k",ls=":")
        axarr2[i,1].axvspan(dr_rsp3ds[0]-dr_rsp3ds[2],dr_rsp3ds[0]+dr_rsp3ds[1],color=assign_color(i,0),alpha=0.3)

	axarr2[i,2].axvline(x=r200m,c="k",ls=":")
        axarr2[i,2].axvspan(dr_rsp3ds[0]-dr_rsp3ds[2],dr_rsp3ds[0]+dr_rsp3ds[1],color=assign_color(i,0),alpha=0.3)
  
	#Plot 3D derivative MCMC fits
	axarr2[i,1].plot(rrange,dr_rhodevs[:,0],linestyle='-',c=assign_color(i,0))
        axarr2[i,1].fill_between(rrange,dr_rhodevs[:,0]-dr_rhodevs[:,2],dr_rhodevs[:,0]+dr_rhodevs[:,1],facecolor=assign_color(i,0),alpha=0.3)
	
        #Plot 3D 2nd derivative MCMC fits
	axarr2[i,2].plot(rrange,dr_rhodevs2[:,0],linestyle='-',c=assign_color(i,0))
        axarr2[i,2].fill_between(rrange,dr_rhodevs2[:,0]-dr_rhodevs2[:,2],dr_rhodevs2[:,0]+dr_rhodevs2[:,1],facecolor=assign_color(i,0),alpha=0.3)

	i+=1


    #Get plots into shape
    axarr1[0,0].set_xscale("log")
    axarr1[1,0].set_xscale("log")
    axarr1[2,0].set_xscale("log")
    axarr1[0,1].set_xscale("log")
    axarr1[1,1].set_xscale("log")
    axarr1[2,1].set_xscale("log")
    axarr2[0,0].set_xscale("log")
    axarr2[1,0].set_xscale("log")
    axarr2[2,0].set_xscale("log")
    axarr2[0,1].set_xscale("log")
    axarr2[1,1].set_xscale("log")
    axarr2[2,1].set_xscale("log")
    axarr2[0,2].set_xscale("log")
    axarr2[1,2].set_xscale("log")
    axarr2[2,2].set_xscale("log")

    axarr1[0,0].set_yscale("log")
    axarr1[1,0].set_yscale("log")
    axarr1[2,0].set_yscale("log")
    axarr2[0,0].set_yscale("log")
    axarr2[1,0].set_yscale("log")
    axarr2[2,0].set_yscale("log")

    axarr1[2,0].set_xlabel(r"$R$ (h$^{-1}$Mpc)")
    axarr1[2,1].set_xlabel(r"$R$ (h$^{-1}$Mpc)")
    axarr2[2,0].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr2[2,1].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr2[2,2].set_xlabel(r"$r$ (h$^{-1}$Mpc)")

    axarr1[0,0].set_ylabel(r"$\xi_{\rm 2D}$($R$)")
    axarr1[1,0].set_ylabel(r"$\xi_{\rm 2D}$($R$)")
    axarr1[2,0].set_ylabel(r"$\xi_{\rm 2D}$($R$)")
    axarr1[0,1].set_ylabel(r"d$\log\xi_{\rm 2D}/$d$\log R$")
    axarr1[1,1].set_ylabel(r"d$\log\xi_{\rm 2D}/$d$\log R$")
    axarr1[2,1].set_ylabel(r"d$\log\xi_{\rm 2D}/$d$\log R$")

    axarr2[0,0].set_ylabel(r"$\xi_{\rm 3D}(r)$")
    axarr2[1,0].set_ylabel(r"$\xi_{\rm 3D}(r)$")
    axarr2[2,0].set_ylabel(r"$\xi_{\rm 3D}(r)$")
    axarr2[0,1].set_ylabel(r"d$\log\xi_{\rm 3D}/$d$\log r$")
    axarr2[1,1].set_ylabel(r"d$\log\xi_{\rm 3D}/$d$\log r$")
    axarr2[2,1].set_ylabel(r"d$\log\xi_{\rm 3D}/$d$\log r$")
    axarr2[0,2].set_ylabel(r"d$^2\log\xi_{\rm 3D}/$d$(\log r)^2$")
    axarr2[1,2].set_ylabel(r"d$^2\log\xi_{\rm 3D}/$d$(\log r)^2$")
    axarr2[2,2].set_ylabel(r"d$^2\log\xi_{\rm 3D}/$d$(\log r)^2$")


    axarr1[0,0].text(0.05, 0.12, r'SNR $42.4$',fontsize=7,transform=axarr1[0,0].transAxes)
    axarr1[1,0].text(0.05, 0.12, r'SNR $43.3$',fontsize=7,transform=axarr1[1,0].transAxes)
    axarr1[2,0].text(0.05, 0.12, r'SNR $30.9$',fontsize=7,transform=axarr1[2,0].transAxes)
    axarr2[0,0].text(0.05, 0.12, r'SNR $42.4$',fontsize=7,transform=axarr2[0,0].transAxes)
    axarr2[1,0].text(0.05, 0.12, r'SNR $43.3$',fontsize=7,transform=axarr2[1,0].transAxes)
    axarr2[2,0].text(0.05, 0.12, r'SNR $30.9$',fontsize=7,transform=axarr2[2,0].transAxes)

    axarr1[0,0].text(0.05, 0.05, r'$M_i-5\log $h$<-19.44$',fontsize=7,transform=axarr1[0,0].transAxes)
    axarr1[1,0].text(0.05, 0.05, r'$M_i-5\log $h$<-18.94$',fontsize=7,transform=axarr1[1,0].transAxes)
    axarr1[2,0].text(0.05, 0.05, r'$M_i-5\log $h$<-18.44$',fontsize=7,transform=axarr1[2,0].transAxes)

    axarr1[0,0].text(2.1,7,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr1[1,0].text(2.1,7,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr1[2,0].text(2.1,7,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr1[0,1].text(2.1,-0.15,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr1[1,1].text(2.1,-0.15,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr1[2,1].text(2.1,-0.15,r'$r_{\mathrm{200m}}$',rotation=90)


    axarr2[0,0].text(0.05, 0.05, r'$M_i-5\log $h$<-19.44$',fontsize=7,transform=axarr2[0,0].transAxes)
    axarr2[1,0].text(0.05, 0.05, r'$M_i-5\log $h$<-18.94$',fontsize=7,transform=axarr2[1,0].transAxes)
    axarr2[2,0].text(0.05, 0.05, r'$M_i-5\log $h$<-18.44$',fontsize=7,transform=axarr2[2,0].transAxes)
    axarr2[0,0].text(2.1,12,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[1,0].text(2.1,12,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[2,0].text(2.1,12,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[0,1].text(2.1,-0.75,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[1,1].text(2.1,-0.75,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[2,1].text(2.1,-0.75,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[0,2].text(2.1,3.2,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[1,2].text(2.1,3.2,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr2[2,2].text(2.1,3.4,r'$r_{\mathrm{200m}}$',rotation=90)

    axarr1[0,0].set_ylim([6e-2,1e1])
    axarr1[1,0].set_ylim([6e-2,1e1])
    axarr1[2,0].set_ylim([6e-2,1e1])
    axarr1[0,1].set_ylim([-2.0,0.0])
    axarr1[1,1].set_ylim([-2.0,0.0])
    axarr1[2,1].set_ylim([-2.0,0.0])
    axarr1[0,0].set_xlim([0,8])
    axarr1[1,0].set_xlim([0,8])
    axarr1[2,0].set_xlim([0,8])
    axarr1[0,1].set_xlim([0,8])
    axarr1[1,1].set_xlim([0,8])
    axarr1[2,1].set_xlim([0,8])

    axarr2[0,0].set_ylim([1e-3,3e1])
    axarr2[1,0].set_ylim([1e-3,3e1])
    axarr2[2,0].set_ylim([1e-3,3e1])
    axarr2[0,1].set_ylim([-3.75,-0.5])
    axarr2[1,1].set_ylim([-3.75,-0.5])
    axarr2[2,1].set_ylim([-3.75,-0.5])
    axarr2[0,2].set_ylim([-3.8,3.8])
    axarr2[1,2].set_ylim([-3.8,3.8])
    axarr2[2,2].set_ylim([-3.8,4])
    axarr2[0,0].set_xlim([0,8])
    axarr2[1,0].set_xlim([0,8])
    axarr2[2,0].set_xlim([0,8])
    axarr2[0,1].set_xlim([0,8])
    axarr2[1,1].set_xlim([0,8])
    axarr2[2,1].set_xlim([0,8])
    axarr2[0,2].set_xlim([0,8])
    axarr2[1,2].set_xlim([0,8])
    axarr2[2,2].set_xlim([0,8])

    axarr1[0,1].xaxis.set_ticklabels([])
    axarr1[0,0].xaxis.set_ticklabels([])
    axarr1[1,0].xaxis.set_ticklabels([])
    axarr1[1,1].xaxis.set_ticklabels([])

    axarr2[0,1].xaxis.set_ticklabels([])
    axarr2[0,0].xaxis.set_ticklabels([])
    axarr2[1,0].xaxis.set_ticklabels([])
    axarr2[1,1].xaxis.set_ticklabels([])
    axarr2[0,2].xaxis.set_ticklabels([])
    axarr2[1,2].xaxis.set_ticklabels([])

    f1.subplots_adjust(hspace=0.1,wspace=0.3)
    f2.subplots_adjust(hspace=0.1,wspace=0.3)

    f1.savefig("%s/plots/2D_graphs_%s.pdf" % (splash_directory, run_type))
    f2.savefig("%s/plots/3D_graphs_%s.pdf" % (splash_directory, run_type))
    plt.close()



def hist_graphs(types):
    f1, axarr1 = plt.subplots(ncols=2,nrows=1,figsize=(6,3))
    f2, axarr2 = plt.subplots(ncols=2,nrows=1,figsize=(6,3))

    i=0
    shifts=np.asarray([1,1,1])
    labels=[r'$M_{\rm i}$-$5\log$h$<-19.44$',r'$M_{\rm i}$-$5\log$h$<-18.94$',r'$M_{\rm i}$-$5\log$h$<-18.44$']
    for type_ in types:
        if mc == False:
            type_ += "_no_mc"
        print("Doing %s",type_)
	dr_data_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (est_directory, type_), sep=' ',header=None,error_bad_lines=False,usecols=(rsteps*4,rsteps*4+1,rsteps*4+3,rsteps*4+4))
	dr_data_array=np.asarray(dr_data_array.values)
	dr_rsp2d_array=dr_data_array[:,0]
	dr_rsp2d_array=dr_rsp2d_array.reshape((dr_rsp2d_array.size,1))
	dr_rsp3d_array=dr_data_array[:,1]
	dr_rsp3d_array=dr_rsp3d_array.reshape((dr_rsp3d_array.size,1))
	dr_maxrhodev_array=dr_data_array[:,2]
	dr_maxinnerrhodev_array=dr_data_array[:,3]
	dr_maxinnerrhodev_array=dr_maxinnerrhodev_array.reshape((dr_maxinnerrhodev_array.size,1))

        print("Data read")
        print("Binning...")

	#Binning and visualizing
	axarr1[0].hist(dr_rsp2d_array.reshape(dr_rsp2d_array.size),100,range=[1,2],color=assign_color(i,0),label=labels[i],density=True,histtype='step')
	axarr1[1].hist(dr_rsp3d_array.reshape(dr_rsp3d_array.size),100,range=[1.2,3],color=assign_color(i,0),label=labels[i],density=True,histtype='step')

	axarr2[0].hist(dr_maxrhodev_array.reshape(dr_maxrhodev_array.size),100,range=[-4.5,-2.5],color=assign_color(i,0),label=labels[i],density=True,histtype='step')
	axarr2[1].hist(dr_maxinnerrhodev_array.reshape(dr_maxinnerrhodev_array.size),100,range=[-6.0,-2.5],color=assign_color(i,0),label=labels[i],density=True,histtype='step')
	i+=1

    print("Plotting...")
    #Putting histograms into shape
    axarr1[0].legend(loc='upper right',fontsize=5)
    axarr2[0].legend(loc='upper right',fontsize=5)

    axarr1[0].set_xlabel(r"$R$ (h$^{-1}$Mpc)")
    axarr1[1].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr2[0].set_xlabel(r"$\frac{\mathrm{d}\log\rho}{\mathrm{d}\log r}(r_{\mathrm{sp}}^{\mathrm{3D}})$")
    axarr2[1].set_xlabel(r"$\frac{\mathrm{d}\log\rho_{\mathrm{in}}f_{\mathrm{trans}}}{\mathrm{d}\log r} (r_{\mathrm{sp}}^{\mathrm{3D}})$")

    axarr1[0].text(2.75, 3.5, r'$R_{\mathrm{sp}}^{\mathrm{2D}}$',fontsize=7)
    axarr1[1].text(2.75, 1.1, r'$r_{\mathrm{sp}}^{\mathrm{3D}}$',fontsize=7)
    axarr1[0].axvline(x=r200m,c="k",linestyle=":")
    axarr1[1].axvline(x=r200m,c="k",linestyle=":")
    axarr1[0].text(1.9,1.5,r'$r_{\mathrm{200m}}$',rotation=90)
    axarr1[1].text(1.9,0.55,r'$r_{\mathrm{200m}}$',rotation=90)
    
    axarr1[0].set_xlim([1.0,3.0])
    axarr1[1].set_xlim([1.0,3.0])

    axarr2[0].set_xlim([-2.5,-5.5])
    axarr2[1].set_xlim([-2.5,-5.5])
    f1.tight_layout()
    f2.tight_layout()
    #plt.tight_layout(w_pad=10)
    f1.savefig("%s/plots/splashback_%s.pdf" % (splash_directory, run_type))
    f2.savefig("%s/plots/derivatives_%s.pdf" % (splash_directory, run_type))
    plt.close()


def comp_splash(types):
    f1 = plt.figure()
    axarr1 = [f1.add_subplot(2, 2, 1), f1.add_subplot(2, 2, 2)]

    types=types
    labels=[r'$M_{\rm i}$-$5\log$h$<-19.44$',r'$M_{\rm i}$-$5\log$h$<-18.94$',r'$M_{\rm i}$-$5\log$h$<-18.44$']
    plt.figure()
    i=0
    yticks=np.linspace(0,0.3,len(types))
    splash2d_mean=np.zeros(0)
    splash2d_low=np.zeros(0)
    splash2d_high=np.zeros(0)
    splash3d_mean=np.zeros(0)
    splash3d_low=np.zeros(0)
    splash3d_high=np.zeros(0)
    colors=np.zeros(0)
    for type_ in types:
        if mc == False:
            type_ += "_no_mc"
        dr_dat=np.genfromtxt("%s/%s/results.csv" % (est_directory, type_), delimiter=',',skip_header=1)
	col=assign_color(i,0)
	lab=labels[i]
        rsp2ds=dr_dat[:,11]
        axarr1[0].errorbar(x=rsp2ds[0],xerr=[[rsp2ds[2]],[rsp2ds[1]]],y=yticks[i],color=col,fmt='s',capsize=4)
        rsp3ds=dr_dat[:,12]
        axarr1[1].errorbar(x=rsp3ds[0],xerr=[[rsp3ds[2]],[rsp3ds[1]]],y=yticks[i]+0.05,color=col,fmt='.',label=lab,capsize=4)

        i+=1       
    f1.legend(fontsize=7,markerscale=1)
    axarr1[0].yaxis.set_ticklabels([])
    axarr1[1].yaxis.set_ticklabels([])
    axarr1[0].yaxis.set_ticks([])
    axarr1[1].yaxis.set_ticks([])
    axarr1[0].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr1[1].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    f1.savefig("%s/plots/splash_comp_%s.pdf" % (splash_directory, run_type))


def color_plot(red_type, blue_type):
    f1, axarr1 = plt.subplots(ncols=2,nrows=1,figsize=(6,3))
    f2, axarr2 = plt.subplots(ncols=3,nrows=1,figsize=(11,3))

    dr_red=np.genfromtxt("%s/%s/xi_2d.dat" % (splash_directory, red_type) )
    dr_blue=np.genfromtxt("%s/%s/xi_2d.dat" % (splash_directory, blue_type) )

    if mc == False:
        red_type += "_no_mc"
        blue_type += "_no_mc"

    dr_full=np.genfromtxt("%s/%s/results.csv" % (est_directory, full_type), delimiter=',',skip_header=1)
    full_rsp2d=dr_full[:,12]
    full_rsp3d=dr_full[:,13]


    idx = (dr_red[:,0] > 0.1) & (dr_red[:,0] < 10.0)
    rr=dr_red[idx,0]
    rrange=np.linspace(rr[0],rr[-1],rsteps)

    red_plotdat=pd.read_csv("%s/%s/plotsave.dat" % (est_directory, red_type),sep=' ',header=None)
    red_plotdat=np.asarray(red_plotdat.values)
    red_sigmas=red_plotdat[0:rsteps,:]
    red_rhos=red_plotdat[rsteps:rsteps*2,:]
    red_devs=red_plotdat[rsteps*2:rsteps*3,:]
    red_rhodevs=red_plotdat[rsteps*3:rsteps*4,:]
    red_rhodevs2=red_plotdat[rsteps*4:rsteps*5,:]


    blue_plotdat=pd.read_csv("%s/%s/plotsave.dat" % (est_directory, blue_type),sep=' ',header=None)
    blue_plotdat=np.asarray(blue_plotdat.values)
    blue_sigmas=blue_plotdat[0:rsteps,:]
    blue_rhos=blue_plotdat[rsteps:rsteps*2,:]
    blue_devs=blue_plotdat[rsteps*2:rsteps*3,:]
    blue_rhodevs=blue_plotdat[rsteps*3:rsteps*4,:]
    blue_rhodevs2=blue_plotdat[rsteps*4:rsteps*5,:]

    red_splashdat=pd.read_csv("%s/%s/values.dat" % (est_directory, red_type), sep=' ',header=None)
    red_splashdat=np.asarray(red_splashdat.values)
    red_rsp2ds=red_splashdat[0,:]
    red_rsp3ds=red_splashdat[1,:]

    blue_splashdat=pd.read_csv("%s/%s/values.dat" % (est_directory, blue_type), sep=' ',header=None)
    blue_splashdat=np.asarray(blue_splashdat.values)
    blue_rsp2ds=blue_splashdat[0,:]
    blue_rsp3ds=blue_splashdat[1,:]

    #Plot 2D Data
    axarr1[0].errorbar(dr_red[:,0], dr_red[:,1], dr_red[:,2], fmt="r.",capsize=2,label =r'$\xi_{DP}^{red}$')
    axarr1[0].errorbar(dr_blue[:,0], dr_blue[:,1], dr_blue[:,2], fmt="b.",capsize=2,label =r'$\xi_{DP}^{blue}$')

    #Plot 2D MCMC fits
    axarr1[0].plot(rrange,red_sigmas[:,0],linestyle='-',c='r')
    axarr1[0].plot(rrange,blue_sigmas[:,0],linestyle='-',c='b')

    axarr1[0].fill_between(rrange,red_sigmas[:,0]-red_sigmas[:,2],red_sigmas[:,0]+red_sigmas[:,1],facecolor='r',alpha=0.3)
    axarr1[0].fill_between(rrange,blue_sigmas[:,0]-blue_sigmas[:,2],blue_sigmas[:,0]+blue_sigmas[:,1],facecolor='b',alpha=0.3)

    #Plot 2D derivative MCMC fits
    axarr1[1].plot(rrange,red_devs[:,0],linestyle='-',c='r')
    axarr1[1].plot(rrange,blue_devs[:,0],linestyle='-',c='b')
    axarr1[1].fill_between(rrange,red_devs[:,0]-red_devs[:,2],red_devs[:,0]+red_devs[:,1],facecolor='r',alpha=0.3)
    axarr1[1].fill_between(rrange,blue_devs[:,0]-blue_devs[:,2],blue_devs[:,0]+blue_devs[:,1],facecolor='b',alpha=0.3)

    #Plot 3D MCMC fits
    axarr2[0].plot(rrange,red_rhos[:,0],linestyle='-',c='r')
    axarr2[0].plot(rrange,blue_rhos[:,0],linestyle='-',c='b')
    axarr2[0].fill_between(rrange,red_rhos[:,0]-red_rhos[:,2],red_rhos[:,0]+red_rhos[:,1],facecolor='r',alpha=0.3)
    axarr2[0].fill_between(rrange,blue_rhos[:,0]-blue_rhos[:,2],blue_rhos[:,0]+blue_rhos[:,1],facecolor='b',alpha=0.3)

    #Plot 3D derivative MCMC fits
    axarr2[1].plot(rrange,red_rhodevs[:,0],linestyle='-',c='r')
    axarr2[1].plot(rrange,blue_rhodevs[:,0],linestyle='-',c='b')

    axarr2[1].fill_between(rrange,red_rhodevs[:,0]-red_rhodevs[:,2],red_rhodevs[:,0]+red_rhodevs[:,1],facecolor='r',alpha=0.3)
    axarr2[1].fill_between(rrange,blue_rhodevs[:,0]-blue_rhodevs[:,2],blue_rhodevs[:,0]+blue_rhodevs[:,1],facecolor='b',alpha=0.3)

    #Plot 3D 2nd derivative MCMC fits
    axarr2[2].plot(rrange,red_rhodevs2[:,0],linestyle='-',c='r')
    axarr2[2].plot(rrange,blue_rhodevs2[:,0],linestyle='-',c='b')

    axarr2[2].fill_between(rrange,red_rhodevs2[:,0]-red_rhodevs2[:,2],red_rhodevs2[:,0]+red_rhodevs2[:,1],facecolor='r',alpha=0.3)
    axarr2[2].fill_between(rrange,blue_rhodevs2[:,0]-blue_rhodevs2[:,2],blue_rhodevs2[:,0]+blue_rhodevs2[:,1],facecolor='b',alpha=0.3)

    #Make vertical line at locations of R200m and R3D
    #axarr[0,0].axvline(x=r200m,c="k",linestyle=":")
    axarr1[0].axvspan(red_rsp2ds[0]-red_rsp2ds[2],red_rsp2ds[0]+red_rsp2ds[1],color='r',alpha=0.3)
    axarr1[0].axvspan(blue_rsp2ds[0]-blue_rsp2ds[2],blue_rsp2ds[0]+blue_rsp2ds[1],color='b',alpha=0.3)
    axarr1[0].axvline(full_rsp2d[0]-full_rsp2d[2],c='k',ls='--')
    axarr1[0].axvline(full_rsp2d[0]+full_rsp2d[1],c='k',ls='--')

    #axarr[0,1].axvline(x=r200m,c="k",linestyle=":")
    axarr1[1].axvspan(red_rsp2ds[0]-red_rsp2ds[2],red_rsp2ds[0]+red_rsp2ds[1],color='r',alpha=0.3)
    axarr1[1].axvspan(blue_rsp2ds[0]-blue_rsp2ds[2],blue_rsp2ds[0]+blue_rsp2ds[1],color='b',alpha=0.3)
    axarr1[1].axvline(full_rsp2d[0]-full_rsp2d[2],c='k',ls='--')
    axarr1[1].axvline(full_rsp2d[0]+full_rsp2d[1],c='k',ls='--')

    #axarr[1,0].axvline(x=r200m,c="k",linestyle=":")
    axarr2[0].axvspan(red_rsp3ds[0]-red_rsp3ds[2],red_rsp3ds[0]+red_rsp3ds[1],color='r',alpha=0.3)
    axarr2[0].axvspan(blue_rsp3ds[0]-blue_rsp3ds[2],blue_rsp3ds[0]+blue_rsp3ds[1],color='b',alpha=0.3)
    axarr2[0].axvline(full_rsp3d[0]-full_rsp3d[2],c='k',ls='--')
    axarr2[0].axvline(full_rsp3d[0]+full_rsp3d[1],c='k',ls='--')

    #axarr[1,1].axvline(x=r200m,c="k",linestyle=":")
    axarr2[1].axvspan(red_rsp3ds[0]-red_rsp3ds[2],red_rsp3ds[0]+red_rsp3ds[1],color='r',alpha=0.3)
    axarr2[1].axvspan(blue_rsp3ds[0]-blue_rsp3ds[2],blue_rsp3ds[0]+blue_rsp3ds[1],color='b',alpha=0.3)
    axarr2[1].axvline(full_rsp3d[0]-full_rsp3d[2],c='k',ls='--')
    axarr2[1].axvline(full_rsp3d[0]+full_rsp3d[1],c='k',ls='--')

    axarr2[2].axvspan(red_rsp3ds[0]-red_rsp3ds[2],red_rsp3ds[0]+red_rsp3ds[1],color='r',alpha=0.3)
    axarr2[2].axvspan(blue_rsp3ds[0]-blue_rsp3ds[2],blue_rsp3ds[0]+blue_rsp3ds[1],color='b',alpha=0.3)
    axarr2[2].axvline(full_rsp3d[0]-full_rsp3d[2],c='k',ls='--')
    axarr2[2].axvline(full_rsp3d[0]+full_rsp3d[1],c='k',ls='--')

    #Trim plots
    axarr1[0].set_xscale("log")
    axarr1[1].set_xscale("log")
    axarr2[0].set_xscale("log")
    axarr2[1].set_xscale("log")
    axarr2[2].set_xscale("log")

    axarr1[0].set_yscale("log")
    axarr2[0].set_yscale("log")

    axarr1[0].set_xlabel(r"$R$ (h$^{-1}$Mpc)")
    axarr1[1].set_xlabel(r"$R$ (h$^{-1}$Mpc)")
    axarr2[0].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr2[1].set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr2[2].set_xlabel(r"$r$ (h$^{-1}$Mpc)")

    axarr1[0].set_ylabel(r"$\xi_{\rm 2D}$($R$)")
    axarr1[1].set_ylabel(r"d$\log\xi_{\rm 2D}/$d$\log R$")
    axarr2[0].set_ylabel(r"$\xi_{\rm 3D}(r)$")
    axarr2[1].set_ylabel(r"d$\log\xi_{\rm 3D}/$d$\log r$")
    axarr2[2].set_ylabel(r"d$^2\log\xi_{\rm 3D}/$d$(\log r)^2$")

    axarr1[0].text(0.125, 0.11, r'$M_i-5\log $h$<-18.94$',fontsize=7)

    f1.subplots_adjust(hspace=0.25,wspace=0.4)
    f2.subplots_adjust(hspace=0.25,wspace=0.4)

    axarr1[0].set_xlim([0,8])
    axarr1[1].set_xlim([0,8])
    axarr2[0].set_xlim([0,8])
    axarr2[1].set_xlim([0,8])
    axarr2[2].set_xlim([0,8])
    
    axarr1[0].set_ylim([0.08,20])

    f1.savefig("%s/plots/color_separated_2D_%s.pdf" % (splash_directory, run_type))
    f2.savefig("%s/plots/color_separated_3D_%s.pdf" % (splash_directory, run_type))
    plt.close()

    f, axarr = plt.subplots(ncols=1,nrows=1,figsize=(3,3))
    red_fraction=np.genfromtxt("%s/%s/fraction.txt" % (est_directory, red_type))
    blue_fraction=np.genfromtxt("%s/%s/fraction.txt" % (est_directory, blue_type))

    axarr.plot(rrange,red_fraction[:,0],c='r')
    axarr.plot(rrange,blue_fraction[:,0],c='b')

    axarr.fill_between(rrange,red_fraction[:,0]-red_fraction[:,2],red_fraction[:,0]+red_fraction[:,1],facecolor='r',alpha=0.3)
    axarr.fill_between(rrange,blue_fraction[:,0]-blue_fraction[:,2],blue_fraction[:,0]+blue_fraction[:,1],facecolor='b',alpha=0.3)

    axarr.axvline(full_rsp3d[0]-full_rsp3d[2],c='k',ls='--')
    axarr.axvline(full_rsp3d[0]+full_rsp3d[1],c='k',ls='--')

    axarr.set_xscale("log")
    axarr.set_xlabel(r"$r$ (h$^{-1}$Mpc)")
    axarr.set_xlim([0,8])
    axarr.set_ylabel(r"fraction")

    f.savefig("%s/plots/color_fraction_%s.pdf" % (splash_directory, run_type)) 

def get_sigdigits(err):
    if (err==0):
        return 1
    xerr=err*1.0
    i=0
    while(xerr<10):
        i=i+1
        xerr=xerr*10.
    return i

def get_tables(types):
    print("Producing fit Table")
    fp=open("%s/plots/fitting_table_%s.txt" % (splash_directory, run_type),"w")
    sp=open("%s/plots/splashback_table_%s.txt" % (splash_directory, run_type),"w")
    header=r'gal cat & $\log_{10}(\rho_{\mathrm{s}})$ & $\log_{10}(\alpha)$ & $\log_{10}(r_{\mathrm{s}})$ & $\rho_{\mathrm{0}}$ & $s_{\mathrm{e}}$ & $\log_{10}(r_{\mathrm{t}})$ & $\log_{10}(\beta)$ & $\log_{10}(\gamma)$ & $\chi^2/\nu$ \\'+'\n'
    spheader=r'gal cat'
    for type_ in types:
        spheader+=r' & %s' % type_
    spheader+=r' \\'
    spheader+='\n'
    fp.write(header)
    sp.write(spheader)
    fp.write("\hline \n")
    fp.write("\hline \n")
    rsp2ds=np.zeros((1,3))
    rsp3ds=np.zeros((1,3))
    for type_ in types:
        fp.write(type_+" & ")
        if mc == False:
            type_ += "_no_mc"
        dat=pd.read_csv("%s/%s/results.csv" % (est_directory, type_))

        foo2d=dat.rsp2d.values
        foo3d=dat.rsp3d.values
        rsp2ds=np.vstack((rsp2ds,foo2d))
        rsp3ds=np.vstack((rsp3ds,foo3d))

        value=dat.lrhos.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lalpha.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lrs.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lrho0.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lse.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lrt.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lbeta.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lgamma.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.chisquare.values
        fp.write(r'$%.*f$ \\'%( 3, value[0]/3.0)+'\n')

        fp.write(r'\hline'+'\n')
    rsp2ds=rsp2ds[1:,:]
    rsp3ds=rsp3ds[1:,:]
    sp.write(r"$R_{\mathrm{sp}}^{\mathrm{2D}}$")
    for i in range(len(rsp2ds)):
        value=rsp2ds[i]
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        sp.write(r' & $%.*f_{-%.*f}^{+%.*f}$'%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))
    sp.write(r' \\'+'\n')
    sp.write(r'\hline'+'\n')    
    sp.write(r"$r_{\mathrm{sp}}^{\mathrm{3D}}$")
    for i in range(len(rsp3ds)):
        value=rsp3ds[i]
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        sp.write(r' & $%.*f_{-%.*f}^{+%.*f}$'%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))
    sp.write(r' \\'+'\n')
    sp.write(r'\hline')
    fp.close()
    sp.close

def get_color_table(red_type, blue_type):
    print("Producing fit Table")
    fp=open("%s/plots/color_fitting_table_%s.txt" % (splash_directory, run_type),"w")
    header=r'gal cat & $\log_{10}(\rho_{\mathrm{s}})$ & $\log_{10}(\alpha)$ & $\log_{10}(r_{\mathrm{s}})$ & $\rho_{\mathrm{0}}$ & $s_{\mathrm{e}}$ & $\log_{10}(r_{\mathrm{t}})$ & $\log_{10}(\beta)$ & $\log_{10}(\gamma)$ & $R_{\mathrm{sp}}^{\mathrm{2D}}$ & $r_{\mathrm{sp}}^{\mathrm{3D}}$ & $\chi^2/\nu$ \\'+'\n'
    spheader=r'gal cat'

    types = [red_type, blue_type]

    for type_ in types:
        spheader+=r' & %s' % type_
    spheader+=r' \\'
    spheader+='\n'
    fp.write(header)
    fp.write("\hline \n")
    fp.write("\hline \n")
    
    if mc == False:
        red_type += "_no_mc"
        blue_type += "_no_mc"

    types = [red_type, blue_type]
    for type_ in types:
        fp.write(type_+" & ")
        dat=pd.read_csv("%s/%s/results.csv" % (est_directory, type_))

        value=dat.lrhos.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lalpha.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lrs.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lrho0.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lse.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lrt.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lbeta.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.lgamma.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.rsp2d.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.rsp3d.values
        sig_digitsm=get_sigdigits(value[1])
        sig_digitsp=get_sigdigits(value[2])
        sigdig = max(sig_digitsp, sig_digitsm)
        fp.write(r'$%.*f_{-%.*f}^{+%.*f}$ & '%( sigdig, value[0], sigdig, value[1], sigdig, value[2]))

        value=dat.chisquare.values
        fp.write(r'$%.*f$ \\'%( 3, value[0]/3.0)+'\n')

        fp.write(r'\hline'+'\n')
    fp.close()
#######################################
#Main
######################################


run_type = "with_hard_spline_no_mc"


mc = False
font={'size':18}
font_small={'size':11}
colors=np.linspace(0,1,6)
add="_best_mc"
rsteps=25
r200m=1.825186

splash_directory = "/work/dominik.zuercher/Output/splashpipe"
est_directory = "/work/dominik.zuercher/Output/Mest"

types=['Planck_PS_21','Planck_PS_21.5','Planck_PS_22']
red_type = 'Planck_PS_21.5_red_hard_spline'
blue_type = 'Planck_PS_21.5_blue_hard_spline'
full_type =  'Planck_PS_21.5'





print("Plotting color separated plot")
color_plot(red_type, blue_type)
get_color_table(red_type, blue_type)
print("Plotting curves")
curve_graphs(types)
print("Plotting splashback comparison")
comp_splash(types)
print("Plotting histograms")
#hist_graphs(types)
print("Producing splashback and fitting parameter tables")
get_tables(types)

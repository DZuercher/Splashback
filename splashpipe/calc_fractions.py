import pandas
import numpy as np
import pylab as pl


def get_sigma(dire, colors,Njack):
    sigma = {}
    xi_full={} 
    xi_err_full={}
    ratio={}
    xi_ratio_full={} 
    xi_err_ratio_full={}
    for col in colors:
        diree = dire+col
        df = pandas.read_csv(diree+"/Clu=Norm_Gal=Norm/pairs.dat", delim_whitespace=1, header=None, names=(["Rp", "DD", "AVR", "Area"]))
        Ncluster = np.loadtxt(diree+"/Clu=Norm_Gal=Norm/Ncluster.dat")
        Ncluster = np.repeat(Ncluster, np.unique(df.Rp.values).size)

        df_r = pandas.read_csv(diree+"/Clu=Ran_Gal=Norm/pairs.dat", delim_whitespace=1, header=None, names=(["Rp", "DR", "AVR", "Area"]))
        Ncluster_r = np.loadtxt(diree+"/Clu=Ran_Gal=Norm/Ncluster.dat")
        Ncluster_r = np.repeat(Ncluster_r, np.unique(df_r.Rp.values).size)

        sigma[col] = df.DD.values/df.Area.values/Ncluster - df_r.DR.values/df_r.Area.values/Ncluster_r


        Nrad=np.unique(df_r.Rp.values).size      
        xi=np.zeros(np.unique(Nrad))
        xi_err = np.zeros(Nrad)
        xi_cov = np.zeros(Nrad*Nrad).reshape(Nrad, Nrad)

	for i in range(Nrad):
	    xi[i] = np.mean(sigma[col][Njack*i:Njack*(i+1)])

	for i in range(Nrad):
	    begi = Njack*i
	    endi = Njack*(i+1)
	    for j in range(Nrad):
		begj = Njack*j
		endj = Njack*(j+1)
		xi_cov[i][j] = np.mean((sigma[col][begi:endi]-xi[i])*(sigma[col][begj:endj]-xi[j]))*(Njack-1.)
	    xi_err[i] = xi_cov[i][i]**0.5
        #np.savetxt("xi_2d%s.dat"% col, np.transpose([10.0**np.unique(df_r.Rp.values), xi, xi_err]))
        #np.savetxt("xi_2d%s_cov.dat"%col, xi_cov)
        xi_full[col]=xi
        xi_err_full[col]=xi_err


    ratio["_red"]=sigma["_red"]/(sigma["_red"]+sigma["_blue"])
    ratio["_blue"]=sigma["_blue"]/(sigma["_red"]+sigma["_blue"])
    for col in colors:
        diree = dire+col
	xi_ratio=np.zeros(np.unique(Nrad))
	xi_err_ratio = np.zeros(Nrad)
	xi_cov_ratio = np.zeros(Nrad*Nrad).reshape(Nrad, Nrad)
	for i in range(Nrad):
	    xi_ratio[i] = np.mean(ratio[col][Njack*i:Njack*(i+1)])

	for i in range(Nrad):
	    begi = Njack*i
	    endi = Njack*(i+1)
	    for j in range(Nrad):
		begj = Njack*j
		endj = Njack*(j+1)
		xi_cov_ratio[i][j] = np.mean((ratio[col][begi:endi]-xi_ratio[i])*(ratio[col][begj:endj]-xi_ratio[j]))*(Njack-1.)
	    xi_err_ratio[i] = xi_cov_ratio[i][i]**0.5
	np.savetxt(diree+"xi_2d_ratio.dat", np.transpose([10.0**np.unique(df_r.Rp.values), xi_ratio, xi_err_ratio]))
	np.savetxt(diree+"xi_2d_ratio_cov.dat", xi_cov_ratio)
        xi_ratio_full[col]=xi_ratio
        xi_err_ratio_full[col]=xi_err_ratio

    return 10.0**np.unique(df_r.Rp.values), xi_full, xi_err_full, xi_ratio_full, xi_err_ratio_full


if __name__ == "__main__":
    directory = "/work/dominik.zuercher/Output/splashpipe/Planck_PS_21.5_blue_orig"
    colors = ["_red", "_blue"]
    Njack=30

    rr, sigma, err, ratio, ratio_err = get_sigma(directory, colors,Njack)

    ax = pl.subplot(221)
    ax.set_xscale("log")
    ax.errorbar(rr, sigma["_red"], yerr=err['_red'], fmt=".r",capsize=2)
    ax.errorbar(rr, sigma["_blue"], yerr=err['_blue'], fmt=".b",capsize=2)
    #ax.plot(rr, sigma["_red"]+sigma["_blue"], color="k")
    ax.set_xlabel(r"$R_{\rm p}$ ($h^{-1}$Mpc)")
    ax.set_ylabel(r"$\Sigma$ ($h^{2}$Mpc$^{-2}$)")
    ax.set_yscale("log")

    ax = pl.subplot(222)
    ax.set_xscale("log")
    ax.errorbar(rr, ratio["_red"], yerr=ratio_err["_red"], fmt=".r",capsize=2)
    ax.errorbar(rr, ratio["_blue"], yerr=ratio_err["_blue"], fmt=".b",capsize=2)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$R_{\rm p}$ ($h^{-1}$Mpc)")
    ax.set_ylabel(r"Fraction")

    pl.tight_layout()
    pl.savefig(directory+"ratio_plot.pdf")

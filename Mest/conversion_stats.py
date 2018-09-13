#! coding=utf8
import numpy as np
import pandas as pd
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
from py_coda import read_pycoda,mcmc
import argparse


def calc_stats(type_, add, modded = False, nwalkers = 28):
    if modded==False:
        chain = np.loadtxt("%s/chains/%s%s/chainstate_%s%s.txt" % (output_dir, type_, add, type_, add))
    else:
        chain = np.loadtxt("%s/chains/%s%s/chainstate_%s%s_modded.txt" % (output_dir, type_, add, type_, add))
    chain = chain.reshape((steps, nwalkers, ndim))
    flatchain = chain.reshape((-1, ndim))
    np.savetxt("%s/chains/%s%s/flatchain_%s%s.txt" % (output_dir, type_, add, type_, add), flatchain)

    params = np.loadtxt("paramnames.in", dtype = str)
    print("Read")
    for num in range(nwalkers):
	chain = flatchain[num::20,:]
	mcmcobj = mcmc(params, chain, thin = 1)

	mcmcobj.geweke(backend = "Agg", savefig = "%s/output/%s%s/%s%s_%s_" % (output_dir, type_, add, type_, add, num),fp="%s/output/%s%s/%s%s_%s_geweke.dat" % (output_dir, type_, add, type_, add, num))

	statfile = open("%s/output/%s%s/%s%s_%s_stats.dat" % (output_dir, type_, add, type_, add, num), "w+")
	mcmcobj.get_stats(fp = statfile)
	statfile.close()

	mcmcobj.plot_traces(backend = "Agg", savefig = "%s/output/%s%s/%s%s_%s_" % (output_dir, type_, add, type_, add, num))
	mcmcobj.heidelberger_welch(fp = "%s/output/%s%s/%s%s_%s_heidelberg_welch.dat" % (output_dir, type_, add, type_, add, num))
	mcmcobj.plot_autocorr(backend = "Agg", savefig = "%s/output/%s%s/%s%s_%s_" % (output_dir, type_, add, type_, add))



if __name__ == "__main__":

    output_dir = "/work/dominik.zuercher/Output/Mest"

    steps = 100000
    ndim = 8
    nwalkers  = 28

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--type_", help="Type",default="Planck_PS_21")
    parser.add_argument("--add", help="Prior",default="")
    parser.add_argument("--modded", help="modded or not",default=0)

    args = parser.parse_args()

    calc_stats(args.type_,args.add, args.modded == 1)


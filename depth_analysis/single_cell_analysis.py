#! coding=utf8
from scipy.stats import linregress
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import sys
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools


def lin_func(x, coeff):
    return np.add(np.multiply(x, coeff[1]), coeff[0])


if __name__ == "__main__": 
    single_cell_path = "/work/dominik.zuercher/singlecell"
    outpath = "./graphs"

    num = sys.argv[1]

    table = np.genfromtxt("%s/cell_%s" % (single_cell_path, num), dtype = str, delimiter = ',')[1:,:]

    #Through away undefined measurments
    select3 = np.where(table[:,7].astype(float) != -999.)[0]
    table = table[select3,:]
    select4 = np.where(table[:,14].astype(float) !=- 999.)[0]
    table = table[select4,:]

    print("The median of the magnitude is: %s" % +str(np.median(table[:,7].astype(np.float))))

    #Through away multiply counted objects
    pastid = 0
    del_pos = np.zeros(0)
    for i in range(0, table.shape[0]):
        objid = table[i,0]
        if objid == pastid:
            del_pos = np.append(del_pos, i)
        pastid = objid
    table = np.delete(table, (del_pos), axis=0)
    table = table.astype(np.float)

    #Separate into galaxy matches and non-matches
    nonfit_pos = tools.vec_hextest(list = table[:,11], hexvalues = flag1success, mode="neg")
    fit_pos = tools.vec_hextest(list = table[:,11], hexvalues = flag1success, mode = "pos")
    fit_exp = table[fit_pos,14]
    nonfit_exp = table[nonfit_pos,14]
    fit_mag = table[fit_pos,7].astype(float)
    nonfit_mag = table[nonfit_pos,7].astype(float)

    drop_index = np.zeros(0)
    drop_bin = np.zeros(0)
    for bin in range(0, np.max(table[:,14]).astype(int)):
        act_exp_fit = np.zeros(0)
        act_exp_nonfit = np.zeros(0)
        act_mag_fit = np.zeros(0)
        act_mag_nonfit = np.zeros(0)

        foo_pos = np.where((fit_exp < bin + 1) & (fit_exp >= bin))
        act_mag_fit = fit_mag[foo_pos]
        foo_pos = np.where((nonfit_exp < bin + 1) & (nonfit_exp >= bin))
        act_mag_nonfit = nonfit_mag[foo_pos]
        tot, binedges = np.histogram(np.append(act_mag_fit, act_mag_nonfit), bins = 100, range = [0,np.max(table[:,7]).astype(int)])
        
        binsize = np.max(table[:,7]).astype(int)/100.
        tot = np.asarray(tot)
        tot = np.log10(tot)
        bincenters = np.add(binedges, 0.5*binsize)
        bincenters = bincenters[:-1]
        
        if ((bin >= 5) & (bin <= 40)):
            xfitdata = bincenters[((bincenters <= 21) & (bincenters >= 15))]
            idx = np.where(((bincenters <= 21) & (bincenters >= 15)))
            yfitdata = tot[idx]
            slope, intersect, r, p, std = linregress(xfitdata, yfitdata)
            devs = np.subtract(tot, lin_func(bincenters, [intersect, slope - std]))
            index = -99
            for i, dev in enumerate(reversed(devs)):
                if dev >= 0:
                    index = len(bincenters) - i
                    break
            if (index == -99):
                print("Failed to find max for bin %s" % bin)
        
        #Plotting
        figure = plt.figure(bin)
        plt.plot(bincenters, tot, 'k.', label = 'All Objects', alpha = 0.7)
        if ((bin >= 5) & (bin <= 40)):
            xp = np.linspace(15, 25, 100)
            yp = lin_func(xp,[intersect, slope])
            plt.plot(xp, yp, 'b-')
            drop_index = np.append(drop_index, bincenters[index])
            drop_bin = np.append(drop_bin, bin)
            plt.axvline(bincenters[index])
        plt.title("Detected object number with %s exposures" % bin)	
        plt.xlabel("mag")
        plt.ylabel("log10(#)")
        figure.savefig("%s/cell_%s/graph_%s.png" % (outpath, num, bin))


    figure = plt.figure()
    plt.plot(drop_bin, drop_index, 'k.')
    plt.xlabel("exposures")
    plt.ylabel("mag")
    plt.ylim(ymin = 20)
    figure.savefig("analysis_%s.png" % num)	 

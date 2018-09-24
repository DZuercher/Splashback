import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def cut(red, ax = 0.8):
    return ax + ((red - 0.1)*0.8/0.23)

if __name__ == "__main__":

    bad_input = "./bad_output.dat"
    good_input = "/work/dominik.zuercher/Output/match_PS/matched.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS"

    """  
    data=np.loadtxt(bad_input)
    redshift=data[:,0]
    color=data[:,1]
    idx=color<900.
    idy=color>-900.
    idn=np.logical_and(idx,idy)
    bad_color=color[idn]
    bad_redshift=redshift[idn]
    """

    data = np.loadtxt(good_input)
    redshift = data[:,0]
    color = data[:,1]
    idx = color < 900.
    idy = color > -900.
    idn = np.logical_and(idx, idy)
    color = color[idn]
    redshift = redshift[idn]


    #plt.scatter(bad_redshift,bad_color,s=0.01,marker='.',color='b')
    plt.scatter(redshift, color, s = 0.01, marker = '.', color = 'r', label = "matched Red galaxies")
    rrange = np.linspace(0, 0.35, 1000)

    yy0 = cut(rrange,0.65)
    yy1 = cut(rrange,0.8)
    yy2 = cut(rrange,0.7)
    yy3 = cut(rrange,0.6)
    yy4 = cut(rrange,0.5)
    yy5 = cut(rrange,0.4)
    yy6 = cut(rrange,0.3)
    yy7 = cut(rrange,0.2)

    plt.plot(rrange, yy0, 'b-', label="orig. cut")
    plt.plot(rrange, yy1, 'b--', label="cut 1")
    plt.plot(rrange, yy2, 'b--', label="cut 2", alpha=0.9)
    plt.plot(rrange, yy3, 'b--', label="cut 3", alpha=0.8)
    plt.plot(rrange, yy4, 'b--', label="cut 4", alpha=0.7)
    plt.plot(rrange, yy5, 'b--', label="cut 5", alpha=0.6)
    plt.plot(rrange, yy6, 'b--', label="cut 6", alpha=0.5)
    plt.plot(rrange, yy7, 'b--', label="cut 7", alpha=0.4)

    plt.ylim([0,3])
    plt.xlim([0,0.33])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("g-r")
    plt.savefig("%s/plot.pdf" % output_dir)

import matplotlib 
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def cut(red, ax = 0.8):
    return ax + ((red - 0.1)*0.8/0.23)

if __name__ == "__main__":

    catalog = "/work/dominik.zuercher/Output/match_PS/matched_spec_new.dat"
    output_dir = "/work/dominik.zuercher/Output/match_PS"


    data = np.loadtxt(catalog)
    redshift = data[:,0]
    red = data[:,1]
    green = data[:,2]
    red_PS = data[:,3]
    green_PS = data[:,4]
    id_ = data[:,5]
    id_PS = data[:,6]


    rrange = np.linspace(0, 0.35, 1000)
    yy0 = cut(rrange, 0.65)

    #Calculate contaminations for mixed characterization (best estimate?)
    #Falses
    reds_below_cut = (id_ == False) & (green_PS - red_PS < 0.55 + redshift*((redshift - 0.1)*12) )
    blues_above_cut = (id_ == True) & (green_PS - red_PS >= 0.55 + redshift*((redshift - 0.1)*12) )
    #Rights
    reds_above_cut = (id_ == False) & (green_PS - red_PS >= 0.55 + redshift*((redshift - 0.1)*12) )
    blues_below_cut = (id_ == True) & (green_PS - red_PS < 0.55 + redshift*((redshift - 0.1)*12) )

    print("Reds below cut: %s" % np.sum(reds_below_cut))
    print("Blues below cut: %s" % np.sum(blues_below_cut))
    print("Reds above cut: %s" % np.sum(reds_above_cut))
    print("Blues above cut: %s" % np.sum(blues_above_cut))
    plt.figure(1)
    plt.scatter(redshift[id_ == True], green[id_ ==True] - red[id_ ==True], s = 0.01, color = 'b', label = "Blue")
    plt.scatter(redshift[id_ == False], green[id_ == False] - red[id_==False], s = 0.01, color = 'r', label = "Red")

    plt.plot(rrange, yy0, 'k-',lw=0.3, label="orig. cut")

    plt.ylim([0,3])
    plt.xlim([0,0.33])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("g-r")
    plt.savefig("%s/plot_1_SDSS_colors.pdf" % output_dir)

   
    plt.figure(2)

    plt.scatter(red[id_ == True], green[id_ == True] - red[id_ == True], s = 0.01, color = 'b', label = "Blue")
    plt.scatter(red[id_ == False], green[id_ == False] - red[id_ == False], s = 0.01, color = 'r', label = "Red")
    plt.xlabel("r")
    plt.ylabel("g-r")
    plt.xlim([-25,-15])
    plt.ylim([-1,2])
    plt.savefig("%s/plot_2_SDSS_colors.pdf" % output_dir)
    
    plt.figure(3)
    plt.scatter(redshift[id_PS == True], green_PS[id_PS ==True] - red_PS[id_PS ==True], s = 0.01, color = 'b', label = "Blue")
    plt.scatter(redshift[id_PS == False], green_PS[id_PS == False] - red_PS[id_PS==False], s = 0.01, color = 'r', label = "Red")
    rrange = np.linspace(0, 0.35, 1000)

    
    yy0 = cut(rrange,0.65)

    plt.plot(rrange, yy0, 'k-',lw=0.3, label="orig. cut")

    plt.ylim([0,3])
    plt.xlim([0,0.33])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("g-r")
    plt.savefig("%s/plot_1_PS_colors.pdf" % output_dir)
    
   
    plt.figure(4)

    plt.scatter(red_PS[id_PS == True], green_PS[id_PS == True] - red_PS[id_PS == True], s = 0.01, color = 'b', label = "Blue")
    plt.scatter(red_PS[id_PS == False], green_PS[id_PS == False] - red_PS[id_PS == False], s = 0.01, color = 'r', label = "Red")
    plt.xlabel("r")
    plt.ylabel("g-r")
    plt.xlim([14,21])
    plt.ylim([-2.5,2.5])
    plt.savefig("%s/plot_2_PS_colors.pdf" % output_dir)
    

    plt.figure(5)
    plt.scatter(redshift[id_ == True], green_PS[id_ ==True] - red_PS[id_ ==True], s = 0.01, color = 'b', label = "Blue")
    plt.scatter(redshift[id_ == False], green_PS[id_ == False] - red_PS[id_==False], s = 0.01, color = 'r', label = "Red")
    rrange = np.linspace(0, 0.35, 1000)

    
    yy0 = cut(rrange,0.65)

    plt.plot(rrange, yy0, 'k-',lw=0.3, label="orig. cut")

    plt.ylim([0,3])
    plt.xlim([0,0.33])
    plt.legend()
    plt.xlabel("z")
    plt.ylabel("g-r")
    plt.savefig("%s/plot_1_mix.pdf" % output_dir)

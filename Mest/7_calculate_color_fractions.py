import numpy as np
import pandas as pd


def calc_color_frac(red_type, blue_type, add, output_dir):

    if mc == True:
	add += "_mc"
	print("Sort according to miscentering parameters...")
	red_parameter_array=pd.read_csv("%s/%s/chainstate_full_modded.txt" % (output_dir, red_type), sep=' ',header=None,error_bad_lines=False,usecols=[8,9])

	blue_parameter_array=pd.read_csv("%s/%s/chainstate_full_modded.txt" % (output_dir, blue_type), sep=' ',header=None,error_bad_lines=False,usecols=[8,9])


	red_parameter_array=np.asarray(red_parameter_array)
	blue_parameter_array=np.asarray(blue_parameter_array)

	print("Ordering")
	red_order=np.lexsort((red_parameter_array[:,0],red_parameter_array[:,1]))
	blue_order=np.lexsort((blue_parameter_array[:,0],blue_parameter_array[:,1]))

    print("Reading Data arrays...")
    red_data_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (output_dir, red_type), sep=' ',header=None,error_bad_lines=False,usecols=range(rsteps,rsteps*2))

    blue_data_array=pd.read_csv("%s/%s/Data_full_modded.txt" % (output_dir, blue_type), sep=' ',header=None,error_bad_lines=False,usecols=range(rsteps,rsteps*2))

    red_rho_array=np.asarray(red_data_array.values)
    blue_rho_array=np.asarray(blue_data_array.values)

    if mc == True:
	print("Reordering arrays...")
	red_rho_array=red_rho_array[red_order,:]
	blue_rho_array=blue_rho_array[blue_order,:]

    red_frac=red_rho_array[red_cutoff:]/(blue_rho_array[blue_cutoff:]+red_rho_array[red_cutoff:])
    blue_frac=blue_rho_array[blue_cutoff:]/(blue_rho_array[blue_cutoff:]+red_rho_array[red_cutoff:])

    red_fracs = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(red_frac, [16, 50, 84], axis=0))))
    blue_fracs = np.asarray(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(blue_frac, [16, 50, 84], axis=0))))

    np.savetxt("%s/%s/fraction.txt" % (output_dir, red_type), red_fracs)
    np.savetxt("%s/%s/fraction.txt" % (output_dir, blue_type), blue_fracs)



if __name__ == "__main__":


    mc = False
    red_type = "Planck_PS_21.5_red_spline_no_mc"
    blue_type = "Planck_PS_21.5_blue_spline_no_mc"
    add = "_best"
    output_dir = "/work/dominik.zuercher/Output/Mest"


    rsteps = 25
    nwalkers = 28
    blue_cutoff = 0#Due to modding files have different length -> cutoff portion at beginning of blue chanis
    red_cutoff =179#Due to modding files have different length -> cutoff portion at beginning of blue chains

    calc_color_frac(red_type, blue_type, add, output_dir)


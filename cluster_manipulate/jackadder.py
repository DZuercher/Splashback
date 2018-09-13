# Adds a Jackknife bin to each cluster in the catalog. 
# Use for normal and random cluster sample

import pandas
import numpy as np


if __name__ == "__main__":

    Njack = 102 #Number of Jackknife bins
    catalog_in = "catalog.dat"
    catalog_out = "new_catalog.csv"


    data = pandas.read_csv(catalog_in, sep = ' ')
    number = data.RA.values.size
    jacks = np.random.randint(Njack, size = number)
    data["jackreg"] = jacks
    data.to_csv(catalog_out, index = False)

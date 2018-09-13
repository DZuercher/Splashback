#Produces an end-user mask meant for distribution 
import sys
sys.path.insert(0,'/home/dominik.zuercher/Documents/RSP_Pro/toolbox')
import tools
import healpy as hp
import numpy as np
import astropy.wcs as wcs
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import frogress
import argparse


def part2(rank):

    steps = 2643 - 635
    chunky = int(steps/size)
    rest = steps - chunky*size
    mini = chunky*rank
    maxi = chunky*(rank + 1)
    if rank >= (size - 1) - rest:
        maxi += 2 + rank - size + rest)
        mini += rank - size + 1 + rest
    if rank == size - 1:
        maxi = steps + 1
    mini += 635
    maxi += 635
    mini = int(mini)
    maxi = int(maxi)


    bad_pixels = np.zeros(0)
    for pro in frogress.bar(range(mini, maxi)):
        flag = False
        for i in range(0, 100):
             try:
                 foohdu = fits.open("%s/mask_%04d.%3d0.fits" % (mask_directory, pro, i), ignore_missing_end = True)
                 break
             except:
                 pass
             if i == 99:
                 flag = True
        if flag == True:
            continue
        w = wcs.WCS(foohdu[1].header)
        for sky in range(0, 100):
            try:
                mask = fits.open('%s/mask_%04d.%03d.fits' % (mask_directory, pro, sky), ignore_missing_end = True)
                body = mask[1].data
                width = body.shape[1]
                height = body.shape[0]
            except:
                continue
            #cut off borders of skycell
            body = body[20:-20,20:-20]
            badpix = np.where(body != 0)
            bady = badpix[0]
            badx = badpix[1]
            bady += 20
            badx += 20
            world = w.all_pix2world(badx, bady, 0, ra_dec_order = True)
            ra = np.radians(world[0])
            dec = np.radians(np.subtract(90.0, world[1]))
            pix = hp.pixelfunc.ang2pix(NSIDE, dec, ra)
            pix = np.unique(pix)
            pix = pix.astype(int)
            bad_pixels = np.append(bad_pixels, pix)

    bad_pixels = np.unique(bad_pixels)
    np.savetxt("%s/mask_pixel_save_%s.dat" % (mid_directory, rank), bad_pixels)

def part3():
    #Part III: combine everything into full maps
    types = ['21','21.5','22']
    for type_ in types: 
        cutmask = hp.read_map("%s/cutmask_bads_%s.fits" % (input_path, type_))
        mask = np.genfromtxt("%s/bad_mask_pixels.dat" % input_path, dtype = int)
        for idx in mask:
            cutmask[idx] = 0
        hp.write_map("%s/full_mask_%s.fits" % (out_path, type_), cutmask, overwrite = True)

        print("Visualization and producing Healpy map...")
        fig = plt.figure(1)
        hp.mollview(cutmask,title = "PS_mask")
        plt.savefig("%s/vis_full_mask_%s.pdf" % (out_path, type_))

if __name__ == "__main__":

    NSIDE = 2**14
    input_path = "."
    mask_directory = "/work/dominik.zuercher/masks"
    out_path = "."


    if part == 1:
        types = ['21','21.5','22']
        for type_ in types: 
            file_ = "%s/bads_%s.fits" % (input_path, type_)
            bads = hp.read_map(file_)
            new_bads = hp.pixelfunc.ud_grade(bads, nside_out = NSIDE)
            hp.write_map("%s/cutmask_bads_%s.fits" % (input_path, type_), new_bads)
    elif part == 2:
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("--mid_directory", help = "Output directory for second part")
        args = parser.parse_args()
        mid_directory = args.mid_directory
        from mpi4py import MPI
        from glob import glob
        comm = MPI.COMM_WORLD
        rank = comm.rank
        size = comm.size
        part2(rank)
        comm.Barrier()
    elif part == 3:
        part3()
    else:
        print("Part parameter is wrong")

#! coding=utf8
#Surhud's Library for getting the coordinates from clicking on picture
import numpy as np
import matplotlib
import matplotlib.image as mpimg
matplotlib.use("TkAgg")
import pylab as pl
import pandas
from glob import glob
import sys
import astropy.wcs as wcs
import numpy as np
from astropy.coordinates import SkyCoord
import csv

def get_ra_dec(x, y, proj):
    c = SkyCoord.from_pixel(x, y, proj, origin=0)
    return c.ra.degree, c.dec.degree

pixscale=0.25/3600

def construct_wcs_alt(ra0, dec0, width, height):
    """Construct world coordinate system."""
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [0.5 * (width + 1), 0.5 * (height + 1)]
    radiusx = 0.5 * width * pixscale
    cdeltx = np.rad2deg(np.tan(np.deg2rad(radiusx))) / (0.5 * (width - 1))
    radiusy = 0.5 * height * pixscale
    cdelty = np.rad2deg(np.tan(np.deg2rad(radiusy))) / (0.5 * (height - 1))
    w.wcs.cdelt = [cdeltx, cdelty]
    w.wcs.crval = [ra0, dec0]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ['deg', 'deg']
    return w


def construct_wcs(ra0, dec0, radius, width):
    proj = wcs.WCS()
    proj.naxis = [width, width]

    proj.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    proj.wcs.cunit = ['deg', 'deg']
    proj.wcs.crval = [ra0, dec0]

    crpix = 0.5 * (width + 1)
    proj.wcs.crpix = [crpix, crpix]
    proj.tunit = ['prob deg-2']
    
    cdelt = np.rad2deg(np.tan(np.deg2rad(radius))) / (0.5 * (width - 1))
    proj.wcs.cdelt = [-cdelt, cdelt]

    proj.wcs.radesys = 'ICRS'

    return proj

def get_params(fname):
    print(fname.split("&"))
    ra, dec, width, height = fname.split("&")
    ra = float(ra.split("=")[1])
    dec = float(dec.split("=")[1])
    width = int(width.split("=")[1])
    height = int(height.split("=")[1][:-5])
    return ra, dec, width, height

class visualizer:

    def __init__(self, itrn=0, fname="debug.dat"):
        self.itrn = itrn
        self.wcs = None
	self.flag=False
        self.fout = open(fname, "a")
        self.fig, self.ax = pl.subplots()
        self.fnames = glob("./Pan-STARRS_Pics/ready_pictures/*")
        self.fig.canvas.callbacks.connect('button_press_event', self.mousecallback)
        self.fig.canvas.callbacks.connect('key_press_event', self.keycallback)
        self.draw()
        pl.show()
    
    def draw(self):
	self.flag=False
        self.ax.clear()
        data = mpimg.imread(self.fnames[self.itrn],format='jpeg')
        data.setflags(write=1)
	try: 
	    data[:,:,0] = data[::-1,:, 0]
            data[:,:,1] = data[::-1,:, 1]
            data[:,:,2] = data[::-1,:, 2]
	except:
	    self.flag=True
	    return
        cenra, cendec, width, height = get_params(self.fnames[self.itrn])
        np.testing.assert_equal(width, height)
        self.wcs = construct_wcs(cenra, cendec, width/2*pixscale, width)
        self.ax.imshow(data, origin="lower")
        self.ax.set_title(r"%.5f %.5f" % (cenra, cendec))
        self.fig.canvas.draw()
        return 
    
    def keycallback(self, event):
        if event.key == "n":
	    if self.flag==True:
	        print("skipping")
	        self.itrn += 1
	        print("Iterator now at: "+str(self.itrn))
                self.draw()
            msg = "%s %d %.7f %.7f\n" % (self.fnames[self.itrn], self.itrn, self.bcgra, self.bcgdec)
            cenra, cendec, width, height = get_params(self.fnames[self.itrn])
            reppos=np.where(abs(lines[1:,0].astype(float)-cenra)<0.00001)[0]
	    if reppos.size!=1:
		print("not working")
            reppos=np.add(reppos,1)
            print(self.bcgra, self.bcgdec)
            lines[reppos,0]=self.bcgra
            lines[reppos,1]=self.bcgdec
            lines[reppos,lines.shape[1]-1]=width
            lines[reppos,lines.shape[1]-2]=height
            self.fout.write(msg)
            sys.stderr.write(msg)
            self.itrn += 1
   	    print("Iterator now at: "+str(self.itrn)) 
        elif event.key == "p":
            self.itrn -= 1
    
        elif event.key == "o":
            msg = "%s %d %.7f %.7f\n" % (self.fnames[self.itrn], self.itrn, self.bcgra, self.bcgdec)
            cenra, cendec, width, height = get_params(self.fnames[self.itrn])
            reppos = np.where(abs(lines[1:,2].astype(float)-cenra) < 0.000001)
            reppos = np.add(reppos,1)
            self.fout.write(msg)
            sys.stderr.write(msg)
            self.bcgra = -99.0
            self.bcgdec = -99.0
            lines[reppos,2] = self.bcgra
            lines[reppos,3] = self.bcgdec
            lines[reppos,lines.shape[1]-1] = width
            lines[reppos,lines.shape[1]-2] = height
            self.itrn += 1
    
        elif event.key == "r":
            self.draw()
        else:
            print "Key action not defined"
    
        self.draw()
    
    def mousecallback(self, event):
        if event.xdata is not None and event.ydata is not None:
            self.bcgra, self.bcgdec = get_ra_dec(event.xdata, event.ydata, self.wcs)
            sys.stderr.write("Mouse picked event : %.2f %.2f" % (event.xdata, event.ydata))
            self.ax.scatter(event.xdata, event.ydata, marker = ".")
            self.fig.canvas.draw()

if __name__ == "__main__":

    inputfile = 'configfile.csv'
    outputfile = 'new_configfile.csv'

    itrn = 0
    fname = 'debug.dat'
    readfile = csv.reader(open(inputfile))
    lines = [l for l in readfile]
    lines = np.asarray(lines)
    lines = np.hstack((lines,np.zeros((lines.shape[0], 2))))
    lines[0,lines.shape[1]-1] = 'HEIGHT'
    lines[0,lines.shape[1]-2] = 'WIDTH'
    writefile = csv.writer(open(outputfile, 'w+'))
    visualizer(itrn, fname)
    writefile.writerows(lines)


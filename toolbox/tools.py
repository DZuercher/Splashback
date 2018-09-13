#! coding=utf8
import scipy.interpolate as spl
import scipy.optimize as op
import pandas as pd
from astropy.io import fits
import pkg_resources
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import frogress
import subprocess as sub
from PIL import Image
import matplotlib
#import matplotlib.image as mpimg
from scipy.ndimage import imread as mpimg
import astropy.wcs as wcs
import emcee
import os
import scipy.integrate as integrate
from scipy.optimize import fmin
from scipy.optimize import fminbound


overlap=240
pixscale = 0.25 #Pixelscale of the Pan-STARRS Survey in arcsec/pix
cosmo = FlatLambdaCDM(H0 = 70 * u.km / u.s / u.Mpc, Om0 = 0.3) #Cosmology used


#Exclusion criteria as inffered from the flags (not actually applied anywhere)
#Flag 1 exclusions
flag1fail=[2048,4096]
#Flag 3 exclusions
flag3fail=[8192]
#Flag 1 galaxy marker (fits to aperature)
flag1success=[2097152]


resource_package = __name__  # Could be any module/package name
resource_path = 'ps1grid.fits' # Do not use os.path.join(), see below
gridfile = pkg_resources.resource_stream(resource_package, resource_path)
gridlist = fits.open(gridfile,ignore_missing_end=True)


def is_number(s):
    """Checks if a string can be converted to a float."""
    try:
         float(s)
         return True
    except ValueError:
         return False


def pixcal(redshift, scalep=1):
    """Calculates the number of pixels corresponding 
    to a certain physical distance if considered at a 
    certain redshift.
    scalep : physical distance in Mpc
    """
    scale = scalep * u.Mpc
    dist = cosmo.luminosity_distance(redshift)
    theta = Angle(np.arctan(np.divide(scale,np.multiply(2, dist))), u.radian)
    totang = np.multiply(2,theta.arcsecond)
    num_of_pix = np.divide(totang,pixscale)
    return num_of_pix.astype(int)


def get_filenames(pros=[],skys=[],ra=[],dec=[],type="stack",filter='i'):
    """Gets filename lines from the Pan-STARRS Archive.
    pros : list of projectioncells
    skys : list of skycells
    lowerlimit / upperlimit : which part of pros and skys to use
    type : the type of the files that should be gathered
    filter : color filter to be used (i,g,r,z,y,all,any)
    """
    pros=np.asarray(pros)
    pros=pros.reshape(pros.size)
    skys=np.asarray(skys)
    skys=skys.reshape(skys.size)
    ra=np.asarray(ra)
    ra=ra.reshape(ra.size)
    dec=np.asarray(dec)
    dec=dec.reshape(dec.size)
    pros=pros.astype(int)
    skys=skys.astype(int)
	
    data=np.zeros((1,9))
    goods = np.zeros(0)
    
    if (len(pros)==0):
        bars=range(ra.size)
    else:
        bars=range(pros.size)

    for ii in frogress.bar(bars):
        if (len(pros)==0):
            url = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra="+str(ra[ii])+"&dec="+str(dec[ii])+"&type="+str(type)
        else:
            url = "http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?skycell="+str(format(pros[ii],'04'))+"."+str(format(skys[ii],'03'))+"&type="+str(type)
        proc = sub.Popen(["wget","-qO-",str(url),"-w 2","--random-wait"],stdout=sub.PIPE)
        dummy = np.zeros((1,9))
        if filter=="all":
            for line in proc.stdout.readlines():
                row=line.split(' ')
                dummy = np.vstack((dummy,row))
            dummy=dummy[2:]
            if (dummy.size!=9 * 5):
                print("Did not get all fitlers for cell or did not find cell")
                continue
        elif filter=='any':
            for line in proc.stdout.readlines():
                row=line.split(' ')
                dummy = np.vstack((dummy,row))
            dummy=dummy[2:]
            if (dummy.size < 9 * 1):
                print("Did not find cell")
                continue
            if (dummy.size < 9 * 5):
                print("Uncomplete")
                cols=dummy[:,4]
                miss=np.setdiff1d(np.asarray(['r','g','i','y','z']),cols)
                add=np.zeros((miss.size,9))
                add=add.astype(str)
                add[:,4]=miss
                dummy=np.vstack((dummy,add))
        else:
            for line in proc.stdout.readlines():
                row=line.split(' ')
                if row[4]==str(filter):
                    dummy = np.vstack((dummy,row))
            dummy=dummy[2:]
            if (dummy.size!=9 * 1):
                print("Did not find cell")
                continue
        data = np.vstack((data, dummy))
        goods = np.append(goods,ii)
    data = data[1:]
    return data,goods.astype(int)


def generate_names(pros=[],skys=[],ra=[], dec=[], redshift=[], scale=1,base="PAN-STARRS",format="jpeg"):
    """Generates a list of names that can be used to name the downloaded pictures from the PAN-STARRS archive which are otherwise very long.
    Have to give ra,dec and redshift or pros,skys for full cells.
    base : First part of the name
    """
    names=np.zeros(0)
    if (len(redshift)!=0):
        redshift=redshift.astype(float)
        wrongi=np.where(redshift==-1.0)
        redshift=np.delete(redshift,wrongi,axis=0)
        ra=np.delete(ra,wrongi,axis=0)
        dec=np.delete(dec,wrongi,axis=0)
        size = pixcal(redshift,scale)
        for i in range(len(ra)):
            name = str(base)+"_ra=" + str(ra[i]) + "&dec=" + str(dec[i]) + "&width=" + str(size[i]) + "&height=" + str(size[i]) + "."+str(format)
            names=np.append(names,name)
    else:
        for i in range(len(pros)):
            name=str(base)+"_"+str(pros[i]).zfill(4)+'.'+str(skys[i]).zfill(3)+"."+str(format)
            names=np.append(names,name)
    return names.reshape(names.size,1)


def generate_https(color1,color2=[],color3=[],pros=[],skys=[],ra=[], dec=[], redshift=[], scale=1,format="jpeg"):
    """Generates a list of URLs that can be used to download the pictures from the PAN-STARRS archive.
    Need to give ra, dec and redshift or pros and skys for the pictures of the full cells.
    format needs to be jpeg or fits
    One color has to be given for fits file and three colors have to be given for jpeg picture. Also if full cells are needed only fits files can be returned.
    Colors have to be the storage paths of the files on the PAN-STARRS archive
    """
    https=np.zeros(0)
    if (len(redshift)==0):
        for i in range(pros.size):
            http="http://ps1images.stsci.edu"+color1[i]
            https=np.append(https,http)
    else:
        redshift=redshift.astype(float)
        wrongi=np.where(redshift==-1.0)
        redshift=np.delete(redshift,wrongi,axis=0)
        ra=np.delete(ra,wrongi,axis=0)
        dec=np.delete(dec,wrongi,axis=0)
        color1=np.delete(color1,wrongi,axis=0)
        color2=np.delete(color2,wrongi,axis=0)
        color3=np.delete(color3,wrongi,axis=0)
        size = pixcal(redshift,scale)
        for i in range(ra.size):
            if format=="jpeg":
                http = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Size=" + str(size[i]) + "," + str(size[i]) + "&Asinh=False&Red=" + str(color1[i]) + "&Green=" + str(color2[i]) + "&Blue=" + str(color3[i]) + "&ra=" + str(ra[i]) + "&dec=" + str(dec[i])
                if (color1[i]=='0')|(color2[i]=='0')|(color3[i]=='0'):
                     id_=np.zeros(3).astype(bool)
                     id_[0] = color1[i]!='0'
                     id_[1] = color2[i]!='0'
                     id_[2] = color3[i]!='0'
                     for j in range(3):
                         if id_[j]==True:
                             if j==0:
                                 http = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Size=" + str(size[i]) + "," + str(size[i]) + "&Asinh=False&Red=" + str(color1[i])+"&ra=" + str(ra[i]) + "&dec=" + str(dec[i])
                             elif j==1:
                                 http = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Size=" + str(size[i]) + "," + str(size[i]) + "&Asinh=False&Red=" + str(color2[i])+"&ra=" + str(ra[i]) + "&dec=" + str(dec[i])
                             else:
                                 http = "http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Size=" + str(size[i]) + "," + str(size[i]) + "&Asinh=False&Red=" + str(color3[i])+"&ra=" + str(ra[i]) + "&dec=" + str(dec[i])
            else:
                http="http://ps1images.stsci.edu"+color1[i]
            https=np.append(https,http)
    return https.reshape(https.size,1)


def get_params(names):
    """Parser extracting ra, dec, width and height from the name of a picture."""
    out=np.zeros((1,4))
    for fname in names:
         ra, dec, width, height = fname.split("&")
         ra = float(ra.split("=")[1])
         dec = float(dec.split("=")[1])
         width = int(width.split("=")[1])
         height = int(height.split("=")[1][:-5])
         out=np.vstack((out,[ra,dec,width,height]))
    out=out[1:,:]
    return out[:,0], out[:,1], out[:,2], out[:,3]


def download(outnames, projectioncells, skycells,ra=[],dec=[],size=0, type='stack',colorcode=['i','r','g'],format_="fits"):
    """Downloads full images (FITS or JPEG) of skycells. If FITS is specified the first color is used
    colorcode : list of three strings being i,g,r,z or y
    format has to be fits or jpeg
    """
    bads=np.zeros(0)
    projectioncells=np.asarray(projectioncells)
    skycells=np.asarray(skycells)
    outnames=np.asarray(outnames)
    for ii in frogress.bar(range(projectioncells.size)):
        if projectioncells.size==1:
            url="http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?skycell="+str(format(int(projectioncells),'04'))+".0"+str(format(int(skycells),'03'))+"&type="+str(type)
        else:
            url="http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?skycell=" + str(projectioncells[ii]) + ".0" + str(skycells[ii])+"&type="+str(type)
        proc = sub.Popen(["wget","-qO-",str(url),"-w 2","--random-wait"],stdout=sub.PIPE)
        colors={"c0" : " ", "c1" : " ", "c2" : " "}
                
        for line in proc.stdout.readlines():
            row=line.split(' ')
            for j in range(len(colorcode)):
                if colorcode[j]==row[4]:
                    colors["c"+str(j)]=row[7]
		
        if colors["c0"]==" ":
            print("Could not locate filenames")
            bads=np.append(bads,ii)
            continue

        if size==0:
            if outnames.size==1:
                if format_=="fits":
                    command = "wget \"http://ps1images.stsci.edu/" + str(colors["c0"]) + "\" -w 2 --random-wait -O "+str(outnames)
                if format_=="jpeg":	
                    command = "wget \"http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Asinh=False&Red=" + str(colors["c0"]) + "&Green=" + str(colors["c1"]) + "&Blue=" + str(colors["c2"]) + "&Size=ALL\" -w 2 --random-wait -O " + str(outnames)
            else:
                if format_=="fits":
                    command = "wget \"http://ps1images.stsci.edu/" + str(colors["c0"]) + "\" -w 2 --random-wait -O "+str(outnames[ii])
                if format_=="jpeg":
                    command = "wget \"http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Asinh=False&Red=" + str(colors["c0"]) + "&Green=" + str(colors["c1"]) + "&Blue=" + str(colors["c2"]) + "&Size=ALL\" -w 2 --random-wait -O " + str(outnames[ii])

        else:
            if outnames.size==1:
                if format_=="fits":
                    command = "wget \"http://ps1images.stsci.edu/" + str(colors["c0"]) + "\" -w 2 --random-wait -O "+str(outnames)
                if format_=="jpeg":
                    command = "wget \"http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Asinh=False&Red=" + str(colors["c0"]) + "&Green=" + str(colors["c1"]) + "&Blue=" + str(colors["c2"]) + "&Size="+str(int(size))+","+str(int(size))+"&ra="+str(ra)+"&dec="+str(dec)+"\" -w 2 --random-wait -O " + str(outnames)
            else:
                if format_=="fits":
                    command = "wget \"http://ps1images.stsci.edu/" + str(colors["c0"]) + "\" -w 2 --random-wait -O "+str(outnames[ii])
                if format_=="jpeg":
                    command = "wget \"http://ps1images.stsci.edu/cgi-bin/fitscut.cgi?Asinh=False&Red=" + str(colors["c0"]) + "&Green=" + str(colors["c1"]) + "&Blue=" + str(colors["c2"]) + "&Size="+str(int(size))+","+str(int(size))+"&ra="+str(ra)+"&dec="+str(dec)+"\" -w 2 --random-wait -O " + str(outnames[ii])

        os.system(command)
    return bads


def edge_check(filenames,infofile):
    """Checks if the pictures whos paths are given as arguments are on the edges of skycells. 
    Returns a dictionary containing the names of the files that are on edge as well as the ra, dec, width, height, projectioncell, sycells 
    and the type of the on edge positioning (outs).
    """
    filenames=np.asarray(filenames)
    ra, dec, width, height = get_params(filenames)
    outs=np.zeros((1,4))
    print("Checking pictures...")
	
    for pic in frogress.bar(range(len(ra))):
        print(filenames[pic])
        data=mpimg(filenames[pic])
        lefts=data[int(height[pic] / 2), 0]
        rights=data[int(height[pic] / 2), int(width[pic])-1]
        bottoms=data[int(height[pic]) - 1, int(width[pic] / 2)]
        tops=data[0, int(width[pic] / 2)]
        leftsum=np.sum(lefts==[255,255,255])
        rightsum=np.sum(rights==[255,255,255])
        bottomsum=np.sum(bottoms==[255,255,255])
        topsum=np.sum(tops==[255,255,255])
        leftout=False
        rightout=False
        bottomout=False
        topout=False
        if leftsum==lefts.shape[0]:
            leftout=True
        if rightsum==rights.shape[0]:
            rightout=True
        if bottomsum==bottoms.shape[0]:	
            bottomout=True
        if topsum==tops.shape[0]:
            topout=True
        outs=np.vstack((outs,[leftout,rightout,bottomout,topout]))

    outs=outs[1:,:]
    edges=(outs[:,0]==True)|(outs[:,1]==True)|(outs[:,2]==True)|(outs[:,3]==True)
    outs=outs[edges]
    filenames=filenames[edges]	
    width=width[edges]
    height=height[edges]
    info = pd.read_csv(infofile)
    pros=np.asarray(info['PRO'].values)
    skys=np.asarray(info['SKY'].values)
    ras=np.asarray(info['RA'].values)
    decs=np.asarray(info['DEC'].values)
    ids=np.zeros(0)

    print("Matching...")
    for id_ in range(ra.size):
        idx=np.where((np.abs(ras-ra[id_])<0.0001) & (np.abs(decs-dec[id_])<0.0001))[0]
        ids=np.append(ids,idx)
    
    ids=ids.astype(int)
    pros=np.take(pros,ids)
    skys=np.take(skys,ids)
    pros=pros[edges]
    skys=skys[edges]
    xcell=np.mod(skys,10)
    ycell=np.divide(skys,10).astype(int)
    ra=ra[edges]
    dec=dec[edges]
    bad1=((xcell==9)&(outs[:,1]==True))
    bad2=((xcell==0)&(outs[:,0]==True))
    bad3=((ycell==0)&(outs[:,2]==True))
    bad4=((ycell==9)&(outs[:,3]==True))
	
    non_fix=((bad1==True)|(bad2==True)|(bad3==True)|(bad4==True))
    fix=np.invert(non_fix)
    print("Search done")
    return filenames[fix], ra[fix], dec[fix], width[fix], height[fix], pros[fix],skys[fix], outs[fix]


def edge_correct(filename,pro,sky,ra,dec,out,width,height,indir,outdir):
    """Complets the given picture and saves it to outname."""
    inname=str(indir)+str(filename)
    outname=str(outdir)+str(filename)
    cenframep=mpimg(inname)

    xcell=int(sky%10)
    ycell=int(sky/10)
    names=np.zeros(0)
    skycells=np.zeros(0)

    if (xcell!=9)&(out[1]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/right.jpeg")
        skycells=np.append(skycells,sky+1)		
    if (xcell!=0)&(out[0]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/left.jpeg")
        skycells=np.append(skycells,sky-1)
    if (ycell!=0)&(out[2]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/bottom.jpeg")
        skycells=np.append(skycells,sky-10)
    if (ycell!=9)&(out[3]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/top.jpeg")
        skycells=np.append(skycells,sky+10)
    if (ycell!=0)&(xcell!=0)&(out[2]==True)&(out[0]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/bottomleft.jpeg")
        skycells=np.append(skycells,sky-11)
    if (ycell!=9)&(xcell!=9)&(out[3]==True)&(out[1]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/topright.jpeg")
        skycells=np.append(skycells,sky+11)
    if (ycell!=0)&(xcell!=9)&(out[2]==True)&(out[1]==True):
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/bottomright.jpeg")
        skycells=np.append(skycells,sky-9)
    if (ycell!=9)&(xcell!=0)&(out[3]==True)&(out[0]==True): 
        names=np.append(names,"./Pan-STARRS_Pics/foo_pics/topleft.jpeg")
        skycells=np.append(skycells,sky+9)

    projectioncells=np.full(names.size,pro)
    if names.size==1:
        names=names[0]
        projectioncells=projectioncells[0]
        skycells=skycells[0]
    bads=download(names,projectioncells,skycells,ra,dec,size=width,format_="jpeg")
    if len(bads)!=0:
        print("skipping")
        return
	
    start,M,zone,xsub,ysub,xsize,ysize,mindec,maxdec,prodec=get_cell_info(pro)
    try:    
        hdulist=fits.open('/work/dominik.zuercher/masks/mask_%04d.%03d.fits' % (pro,sky),ignore_missing_end=True)
    except:
        return
    system=wcs.WCS(hdulist[1].header)
    pix=system.all_world2pix([[ra,dec]],0,ra_dec_order=True)
    x0=pix[:,0]
    y0=ysub-pix[:,1]
	
    x00=np.zeros(0)
    y00=np.zeros(0)
    for i in range(skycells.size):
        if skycells.size==1:
            hdulist=fits.open('/work/dominik.zuercher/masks/mask_%04d.%03d.fits' % (int(projectioncells),int(skycells)),ignore_missing_end=True)
        else:
            hdulist=fits.open('/work/dominik.zuercher/masks/mask_%04d.%03d.fits' % (int(projectioncells[i]),int(skycells[i])),ignore_missing_end=True)
        system=wcs.WCS(hdulist[1].header)
        pix=system.all_world2pix([[ra,dec]],0,ra_dec_order=True)
        x00=np.append(x00,pix[:,0])
        y00=ysub-np.append(y00,pix[:,1])

    while 1==1:
        if (out[0]==1)&(out[2]==1):
            cenframe = cenframep[:-int(height/2-(ysub-y0)+overlap),int(width/2-x0+overlap):]
            try:
                leftframep = mpimg("./Pan-STARRS_Pics/foo_pics/left.jpeg")
            except:
                return
            leftframe = leftframep[:-int(height/2-(ysub-y00[0])+overlap),:-int(width/2-(xsub-x00[0])+overlap)]
            tot = np.hstack((leftframe,cenframe))
            try:
                bottomframep = mpimg("./Pan-STARRS_Pics/foo_pics/bottom.jpeg")
            except:
                return
            bottomframe = bottomframep[int(height/2-y00[1]+overlap):, int(width/2-x00[1]+overlap):]
            try:
                bottomleftframep = mpimg("./Pan-STARRS_Pics/foo_pics/bottomleft.jpeg")
            except:
                return
	    bottomleftframe = bottomleftframep[int(height/2-y00[2]+overlap):,:-int(width/1-(xsub-x00[2])+overlap)]
            bottomrow=np.hstack((bottomleftframe,bottomframe))
            tot=np.vstack((tot,bottomrow))
            break
        if (out[0]==1)&(out[3]==1):
            cenframe = cenframep[int(height/2-y0+overlap):,int(width/2-x0+overlap):]
	    try:
                leftframep = mpimg("./Pan-STARRS_Pics/foo_pics/left.jpeg")
            except:
                return
            leftframe = leftframep[int(height/2-y00[0]+overlap):, :-int(width/2-(xsub-x00[0])+overlap)]
            tot = np.hstack((leftframe,cenframe))
            try:
                topframep = mpimg("./Pan-STARRS_Pics/foo_pics/top.jpeg")
            except:
                return
            topframe = topframep[:-int(height/2-(ysub-y00[1])+overlap),int(width/2-x00[1]+overlap):]
            try:
                topleftframep = mpimg("./Pan-STARRS_Pics/foo_pics/topleft.jpeg")
            except:
                return
	    topleftframe = topleftframep[:-int(height/2-(ysub-y00[2])+overlap),:-int(width/2-(xsub-x00[2])+overlap)]
            toprow=np.hstack((topleftframe,topframe))
            tot=np.vstack((toprow,tot))
            break
        if (out[1]==1)&(out[2]==1):
            cenframe = cenframep[:-int(height/2-(ysub-y0)+overlap),:-int(width/2-(xsub-x0)+overlap)]
            try:
                rightframep = mpimg("./Pan-STARRS_Pics/foo_pics/right.jpeg")
            except:
                return
            rightframe = rightframep[:-int(height/2-(ysub-y00[0])+overlap),int(width/2-x00[0]+overlap):]
            tot=np.hstack((cenframe,rightframe))
            try:
                bottomframep = mpimg("./Pan-STARRS_Pics/foo_pics/bottom.jpeg")
            except:
                return
            bottomframe = bottomframep[int(height/2-y00[1]+overlap):,:-int(width/2-(xsub-x00[1])+overlap)]
	    try:
                bottomrightframep = mpimg("./Pan-STARRS_Pics/foo_pics/bottomright.jpeg")
            except:
                return
            bottomrightframe = bottomrightframep[int(height/2-y00[2]+overlap):, int(width/2-x00[2]+overlap):]
            bottomrow=np.hstack((bottomframe,bottomrightframe))
            tot=np.vstack((tot,bottomrow))
            break
        if (out[1]==1)&(out[3]==1):
            cenframe = cenframep[int(height/2-y0+overlap):,:-int(width/2-(xsub-x0)+overlap)]
            try:
                rightframep = mpimg("./Pan-STARRS_Pics/foo_pics/right.jpeg")
            except:
                return
            rightframe = rightframep[int(height/2-y00[0]+overlap):, int(width/2-(xsub-x00[0])+overlap):]
            tot=np.hstack((cenframe,rightframe))
            try:
                topframep = mpimg("./Pan-STARRS_Pics/foo_pics/top.jpeg")
            except:
                return
            topframe = topframep[:-int(height/2-(ysub-y00[1])+overlap),:-int(width/2-(xsub-x00[1])+overlap)]
	    try:
                toprightframep = mpimg("./Pan-STARRS_Pics/foo_pics/topright.jpeg")
            except:
                return
            toprightframe = toprightframep[:-int(height/2-(ysub-y00[2])+overlap),int(width/2-x00[2]+overlap):]
            toprow=np.hstack((topframe,toprightframe))
            tot=np.vstack((toprow,tot))
            break
        if out[0]==1:
            cenframe = cenframep[:,int(width/2-x0+overlap):]
            try:
                leftframep = mpimg("./Pan-STARRS_Pics/foo_pics/left.jpeg")
            except:
                return
            leftframe = leftframep[:,:-int(width/2-(xsub-x00)+overlap)]
            tot = np.hstack((leftframe,cenframe))
            break
        if out[1]==1:        
            cenframe = cenframep[:,:-int(width/2-(xsub-x0)+overlap)]
            try:
                rightframep = mpimg("./Pan-STARRS_Pics/foo_pics/right.jpeg")
            except:
                return
            rightframe = rightframep[:,int(width/2-x00+overlap):]
            tot=np.hstack((cenframe,rightframe))
            break
        if out[2]==1:
            cenframe = cenframep[:-int(height/2-(ysub-y0)+overlap),:]
            try:
                bottomframep = mpimg("./Pan-STARRS_Pics/foo_pics/bottom.jpeg")
            except:
                return
            bottomframe = bottomframep[int(height/2-y00+overlap):,:]
            tot=np.vstack((cenframe,bottomframe))
            break
        if out[3]==1:
            cenframe = cenframep[int(height/2-y0+overlap):,:]
            try:
                topframep = mpimg("./Pan-STARRS_Pics/foo_pics/top.jpeg")
            except:
                return
	    topframe = topframep[:-int(height/2-(ysub-y00)+overlap),:]
            tot=np.vstack((topframe,cenframe))
            break

    im = Image.fromarray(tot).convert('RGB')
    im.save(str(outname))


def get_cell_info(pros):
    head = gridlist[0].header
    body = gridlist[1].data
    decs=body['DEC']
    idx=np.zeros(0)
    pros=np.asarray(pros)
    if pros.size==1:
        idy=np.where(body['PROJCELL']<=pros)[0]
        idx=np.append(idx,idy[-1])
    else:
        for i in range(pros.size):
            idy=np.where(body['PROJCELL']<=pros[i])[0]
            idx=np.append(idx,idy[-1])	
    idx=idx.astype(int)
    startcell=body['PROJCELL'][idx]
    M=body["NBAND"][idx]
    zone=body['ZONE'][idx]
    xsub=body['XCELL'][idx]
    ysub=body['YCELL'][idx]
    xsize=np.add(np.multiply(np.subtract(xsub,480),10),480)
    ysize=np.add(np.multiply(np.subtract(ysub,480),10),480)
    mindec=body['DEC_MIN'][idx]
    maxdec=body['DEC_MAX'][idx]
    prodec = body['DEC'][idx]
    return startcell,M,zone,xsub,ysub,xsize,ysize,mindec,maxdec, prodec


def get_ra_dec(x, y, proj):
    c = SkyCoord.from_pixel(x, y, proj, origin=0)
    return c.ra.degree, c.dec.degree


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


def get_centers(pros,skys,borders=False,maskdetect=False):
    """If borders is True returns ra,dec borders of the cell instead of the center"""
    starts,Ms,zones,xsubs,ysubs,xsizes,ysizes,mindecs,maxdecs,cdecs=get_cell_info(pros)
	
    minras=np.zeros(0)
    mindecs=np.zeros(0)
    maxras=np.zeros(0)
    maxdecs=np.zeros(0)
    raout=np.zeros(0)
    decout=np.zeros(0)
    pros=np.asarray(pros)
    skys=np.asarray(skys)
    for i in range(pros.size):
        skyx=np.linspace(-xsizes[i]/2.0+xsubs[i]/2.0,xsizes[i]/2.0-xsubs[i]/2.0,10)
        skyy=np.linspace(-ysizes[i]/2.0+ysubs[i]/2.0,ysizes[i]/2.0-ysubs[i]/2.0,10)
        xcen=skyx[np.mod(skys[i],10.0).astype(int)]
        ycen=skyy[np.divide(skys[i],10.0).astype(int)]
        flag=False
        try:
            mask = fits.open('/work/dominik.zuercher/masks/mask_%04d.%03d.fits' % (pros[i],skys[i]),ignore_missing_end=True)
            print("passed")
        except:
            if maskdetect==True:
                print("Mask not found")
                if borders==False:
                    raout=np.append(raout,-99)
                    decout=np.append(decout,-99)
                else:
                    mindecs=np.append(mindecs,-99)
                    minras=np.append(minras,-99)
                    maxdecs=np.append(maxdecs,-99)
                    maxras=np.append(maxras,-99)
                continue

            for j in range(0,100):
                try:
                    mask = fits.open('/work/dominik.zuercher/masks/mask_%04d.%03d.fits' % (pros[i],j),ignore_missing_end=True)
                    break
                except:
                    pass
                if j==99:
                    print("No mask in this projectioncell!")
                    if borders==False:
                        raout=np.append(raout,-99)
                        decout=np.append(decout,-99)
                    else:
                        mindecs=np.append(mindecs,-99)
                        minras=np.append(minras,-99)
                        maxdecs=np.append(maxdecs,-99)
                        maxras=np.append(maxras,-99)
	            flag=True

        if flag==True:
            continue
        try:
            rel1=mask[1].header['CRPIX1']
	    rel2=mask[1].header['CRPIX2']
            w=wcs.WCS(mask[1].header)
        except:
            print("Broken Mask")
            if borders==False:
                raout=np.append(raout,-999)
                decout=np.append(decout,-999)
            else:
                mindecs=np.append(mindecs,-999)
                minras=np.append(minras,-999)
                maxdecs=np.append(maxdecs,-999)
                maxras=np.append(maxras,-999)
            continue
		
        xcen+=rel1
        ycen+=rel2
		
        if borders==True:
            xmin=xcen-xsubs[i]/2.0
            xmax=xcen+xsubs[i]/2.0
            ymin=ycen-ysubs[i]/2.0
            ymax=ycen+ysubs[i]/2.0
        if borders==False:
            foo=w.all_pix2world(xcen,ycen,0,ra_dec_order=True)
            raout=np.append(raout,foo[0])
            decout=np.append(decout,foo[1])
        if borders==True:
            mini=w.all_pix2world(xmin,ymin,0,ra_dec_order=True)
            maxi=w.all_pix2world(xmax,ymax,0,ra_dec_order=True)
            if pros[i]!=2643:
                minras=np.append(minras,mini[0])
                maxras=np.append(maxras,maxi[0])
                mindecs=np.append(mindecs,mini[1])
                maxdecs=np.append(maxdecs,maxi[1])
            else:
                data=get_filenames(pros=np.repeat(int(pros[i]),3),skys=[int(skys[i]),int(skys[i])-1,int(skys[i])+1])
                try:
                    ra1=float(data[0][0][2])
                    ra2=float(data[0][1][2])
                    ra3=float(data[0][2][2])
                    minra=ra3+(ra1-ra3)/2.0
                    maxra=ra1+(ra2-ra1)/2.0
                    minras=np.append(minras,maxra)
                    maxras=np.append(maxras,minra)
                    mindecs=np.append(mindecs,mini[1])
                    maxdecs=np.append(maxdecs,maxi[1])
                except:
                    print("Broken Mask")
                    if borders==False:
                        raout=np.append(raout,-999)
                        decout=np.append(decout,-999)
                    else:
                        mindecs=np.append(mindecs,-999)
                        minras=np.append(minras,-999)
                        maxdecs=np.append(maxdecs,-999)
                        maxras=np.append(maxras,-999)
    if borders==False:
        return raout, decout
    else:
        return maxras,minras,mindecs,maxdecs


def create_grid(pros):
    """Creates a list of projectioncells in first column and all 100 skycells per projectioncell in the second column"""
    fullskys=np.arange(0,100,1)
    fulls=np.zeros((0,2))
    for pro in pros:
        vec1=np.repeat(pro,100)
        vec2=np.hstack((vec1.reshape((100,1)),fullskys.reshape((100,1))))
        fulls=np.vstack((fulls,vec2))
    return fulls


def hextest(element,hexvalues,mode='pos'):
    """Checks if element has a certain hexbit set or not.	
    element and hexvalues are just integer numbers.
    mode can be pos or neg: For pos True is returned if the bit is            set. 
    """
    foo=np.zeros_like(hexvalues,dtype=bool)
    for i,hexvalue in enumerate(hexvalues):
        testbit=int(int(element)/int(hexvalue))
        if (testbit%2!=0):
            foo[i]=True #Meaning that the bit is set
    if mode=='neg':
        if np.any(foo==True):
            return False
        else:
            return True
    if mode=='pos':
        if np.all(foo==True):
            return True
        else:
            return False


def vec_hextest(list,hexvalues,mode):
    """Vectorized version of the hextest function. 
    Takes a list of elements instead of just one integer."""
    outarray=np.zeros_like(list,dtype=bool)
    for i,element in enumerate(list):
        outarray[i]=hextest(element,hexvalues,mode)
    return outarray


def decide_pro(ra,dec):
    dec=dec.reshape((1,dec.size))
    head = gridlist[0].header
    body = gridlist[1].data
        
    mins=body['DEC_MIN']
    maxs=body['DEC_MAX']
    mins=mins.reshape((mins.size,1))
    maxs=maxs.reshape((maxs.size,1))

    idmin= dec<=mins
    idmax= dec>=maxs
    idlowout= np.all(idmin,axis=0)
    idhighout= np.all(idmax,axis=0)

    idtot= idmin==idmax
    where=np.where(idtot==True)
    flat=where[0]
    flatid=where[1]
		
    if np.size(flatid)+np.sum(idlowout)+np.sum(idhighout)!=dec.size:
        print("Error: Something went wrong...")
        print(ra,dec)
	
    zones=np.zeros(dec.size)
    zones[flatid]=flat
    zones[idlowout]=0
    zones[idhighout]=mins.shape[0]-1
    zones=zones.astype(int)
    M=np.take(body['NBAND'],zones)
    idx=np.zeros(0)
    for ii,m in enumerate(M):
        ras=np.asarray([n*360.0/m for n in range(m)])
        idxp=np.argmin(np.abs(np.subtract(ras,ra[ii])))
        if idxp==M[ii]-1:
            if np.abs(np.abs(ra[ii]-ras[0])-360.0) < np.min(np.abs(np.subtract(ras,ra[ii]))) :
                idxp=0 
        idx=np.append(idx,idxp)
    startcell=body['ProjCell'][zones]
    projectioncell=np.add(startcell,idx).astype(int)
    return projectioncell

def decide_sky(ra,dec,cellnumber):
    '''Note that cellnumber is scalar !'''
    tdec=dec[0]
    head = gridlist[0].header
    body = gridlist[1].data

    mins=body['DEC_MIN']
    maxs=body['DEC_MAX']
    mins=mins.reshape((mins.size,1))
    maxs=maxs.reshape((maxs.size,1))

    idmin= tdec<=mins
    idmax= tdec>=maxs
    idlowout= np.all(idmin,axis=0)
    idhighout= np.all(idmax,axis=0)

    idtot= idmin==idmax
    where=np.where(idtot==True)
    flat=where[0]
    flatid=where[1]

    if np.size(flatid)+np.sum(idlowout)+np.sum(idhighout)!=tdec.size:
        print("Error: Something went wrong...")
        print(ra,dec)

    zones=np.zeros(tdec.size)
    zones[flatid]=flat
    zones[idlowout]=0
    zones[idhighout]=mins.shape[0]-1
    zones=zones.astype(int)

    xsub=body['XCELL'][zones]
    ysub=body['YCELL'][zones]
    xsize=(xsub-480)*10+480
    ysize=(ysub-480)*10+480
    #Decide on SkyCell
    skyx=np.linspace(-xsize/2+xsub/2,xsize/2-xsub/2,10)
    skyy=np.linspace(-ysize/2+ysub/2,ysize/2-ysub/2,10)

    names=["/work/dominik.zuercher/masks/mask_"+str(int(cellnumber)).zfill(4)+".0"+str(n).zfill(2)+".fits" for n in range(0,100)]
    flag=False
    for i,name in enumerate(names):
         try:
             foohdu=fits.open(name,ignore_missing_end=True)
             break	
         except:
             pass
         if i==99: #if no masks available for this cell
             print("No mask available in this Projectioncell")
             flag=True

    if flag==True:
        return np.repeat(-99,ra.size)
    try:
        wcsfoo=wcs.WCS(foohdu[1].header)
        rel1=foohdu[1].header['CRPIX1']
        rel2=foohdu[1].header['CRPIX2']
    except:
        return np.zeros(0)
    ra=ra.reshape((ra.size,1))
    dec=dec.reshape((dec.size,1))
    world=np.hstack((ra,dec)).reshape((ra.size,2))
    pix=wcsfoo.all_world2pix(world,0,ra_dec_order=True)
    xpix=np.subtract(pix[:,0],rel1)
    ypix=np.subtract(pix[:,1],rel2)
    if xpix.size==1:
         print("Rare single case")
         xmin=np.argmin(np.abs(skyx-xpix))
         ymin=np.argmin(np.abs(skyy-ypix))
    else:
         fooarrayx=np.array([skyx.reshape(skyx.size,1),]*xpix.size).transpose()
         fooarrayx=fooarrayx[0]
         fooarrayx=np.subtract(fooarrayx,xpix.reshape(1,xpix.size))
         fooarrayx=np.absolute(fooarrayx)
         xmin=np.argmin(fooarrayx,axis=0)
         xmin=xmin.reshape(xmin.size,1)

         fooarrayy=np.array([skyy.reshape(skyy.size,1),]*ypix.size).transpose()
         fooarrayy=fooarrayy[0]
         fooarrayy=np.subtract(fooarrayy,ypix.reshape(1,ypix.size))
         fooarrayy=np.absolute(fooarrayy)
         ymin=np.argmin(fooarrayy,axis=0)
         ymin=ymin.reshape(ymin.size,1)

    skycells=np.add(np.multiply(ymin,10),xmin)
    return skycells


def mask_select(output,cellnumber,id_=[]):
    '''Note that cellnumber is scalar !'''
    usky=np.unique(output[:,3].astype(int))
    selectedobjects=np.zeros((1,output.shape[1]))
    selectedids=np.zeros((1,1))
    for n in range(usky.size):
        skycell=usky[n]
        print("Doing "+str(cellnumber)+"."+str(skycell))
        try:
            mask = fits.open('/work/dominik.zuercher/masks/mask_%04d.%03d.fits' % (cellnumber,skycell),ignore_missing_end=True)
        except:
            print("Mask not found")
            continue
        idy=(output[:,3].astype(int)==skycell)
        validobjects=output[idy,:]
        if len(id_)>=1:
            validid_=id_[idy,:]
        try:
            body=mask[1].data
            width=body.shape[1]
            height=body.shape[0]
            w=wcs.WCS(mask[1].header)
        except:
            continue
        pixval=w.all_world2pix(validobjects[:,0:2],0,ra_dec_order=True)
        off=(pixval[:,1]>=height)|(pixval[:,1]<0)|(pixval[:,0]<0)|(pixval[:,0]>=width)
        if np.sum(off)!=0:
            print("Error: off limit values")
            print(validobjects[off,:])
            on=np.invert(off)
            pixval=pixval[on,:]
            validobjects=validobjects[on,:]
            if len(id_)>=1:
                validid_=validid_[on,:]
        foo=body[pixval[:,1].astype(int),pixval[:,0].astype(int)]
        check=(foo==0)
        selectedobjects=np.vstack((selectedobjects,validobjects[check,:]))
        if len(id_)>=1:  
            selectedids=np.vstack((selectedids,validid_[check,:]))
    if len(id_)>=1:
        selectedobjects=selectedobjects[1:,:]
        selectedids=selectedids[1:,:]
        return selectedobjects,selectedids
    else:
        selectedobjects=selectedobjects[1:,:]
        return selectedobjects



def rho_inner(r,p):
    rhos, alpha, rs = p
    rhos = 10.0**rhos
    alpha = 10.0**alpha
    rs = 10.0**rs
    rhoin = rhos * np.exp(-2.0 / alpha * ((r / rs)**alpha -1.0))
    return rhoin

def rho_outer(r,p):
    rho0, se = p
    #rho0=10.0**rho0
    rhoout = rho0 * (r/1.5)**(-se)
    return rhoout

def ftrans(r,p):
    rt, beta, gamma = p
    rt = 10.0**rt
    beta = 10.0**beta
    gamma = 10.0**gamma
    ftrans = (1.0 + (r / rt)**beta)**(-gamma / beta)
    return ftrans

def einasto_profile(rr, p):
    '''3D galaxy density distribution as proposed by Diemer ] Kravtsov 2014 '''
    p0 = p[0:3] 
    p1 = p[3:5]
    p2 = p[5:8]
    if rr.size==1:
        rho=rho_inner(rr,p0) * ftrans(rr,p2) + rho_outer(rr,p1)
        return rho
    else:
        res=np.zeros(0)
        for r in rr:
            rho=rho_inner(r,p0) * ftrans(r,p2) + rho_outer(r,p1)
            res=np.append(res,rho)
        return res

def einasto_profile_inner(r, p):
    '''3D galaxy density distribution as proposed by Diemer ] Kravtsov 2014 '''
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rho=rho_inner(r,p0) * ftrans(r,p2)
    return rho

def logdev_einasto_profile_inner(r,p):
    rhos, alpha, rs, rho0, se, rt, beta, gamma = p
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rhos=10.0**rhos
    #rho0=10.0**rho0
    alpha = 10.0**alpha
    rs = 10.0**rs
    gamma = 10.0**gamma
    beta = 10.0**beta
    rt = 10.0**rt
    #return r*(ftrans(r,p2)*rho_inner(r,p0))/einasto_profile(r,p)*( -2.0*r**(alpha-1.0)/rs**alpha - gamma*r**(beta-1.0)/rt**beta*(1.0+(r/rt)**beta)**(-1.0)  )
    return 2*(r/rs)**alpha+gamma*1./(1+(r/rt)**beta)*(r/rt)**beta

def einasto_profile_outer(r, p):
    '''3D galaxy density distribution as proposed by Diemer ] Kravtsov 2014 '''
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rho=rho_outer(r,p1)
    return rho

def dev_einasto_profile(r, p):
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rhos, alpha, rs, rho0, se, rt, beta, gamma = p
    rhos=10.0**rhos
    #rho0=10.0**rho0
    alpha = 10.0**alpha
    rs = 10.0**rs
    gamma = 10.0**gamma
    beta = 10.0**beta
    rt = 10.0**rt
    return -rho_inner(r, p0) * ftrans(r, p2) * (2. * r**(alpha-1.)/rs**alpha + gamma * r**(beta-1.)/rt**beta/(1+(r/rt)**beta) ) - se/r * rho_outer(r,p1)

def logdev_einasto_profile(rr, p):
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rhos, alpha, rs, rho0, se, rt, beta, gamma = p
    rhos=10.0**rhos
    #rho0=10.0**rho0
    alpha = 10.0**alpha
    rs = 10.0**rs
    gamma = 10.0**gamma
    beta = 10.0**beta
    rt = 10.0**rt
    result=np.zeros(0)
    if rr.size==1:
        return -rr*rho_inner(rr,p0)*ftrans(rr,p2)/einasto_profile(rr,p) * (2. * rr**(alpha-1.)/rs**alpha + gamma * rr**(beta-1.)/rt**beta*(1+(rr/rt)**beta)**(-1.0) ) - se*rho_outer(rr,p1)/einasto_profile(rr,p)
    for r in rr:
        result=np.append(result,-r*rho_inner(r,p0)*ftrans(r,p2)/einasto_profile(r,p) * (2. * r**(alpha-1.)/rs**alpha + gamma * r**(beta-1.)/rt**beta*(1+(r/rt)**beta)**(-1.0) ) - se*rho_outer(r,p1)/einasto_profile(r,p))
    return result

def dev2_einasto_profile(r, p):
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rhos, alpha, rs, rho0, se, rt, beta, gamma = p
    rhos=10.0**rhos
    #rho0=10.0**rho0
    alpha = 10.0**alpha
    rs = 10.0**rs
    gamma=10.0**gamma
    beta=10.0**beta
    rt=10.0**rt
    return rho_inner(r, p0) * ftrans(r, p2) * ( (2. * r**(alpha - 1.0)/rs**alpha + gamma * r**(beta - 1.0)/rt**beta*(1+(r/rt)**beta)**(-1.0))**2. - 2.*(alpha - 1.) * r**(alpha-2.)/rs**alpha - gamma * (beta - 1.0) * r**(beta - 2.)/rt**beta*(1+(r/rt)**beta)**(-1.0)+gamma*r**(beta-1.0)/rt**beta*beta*r**(beta-1)/rt**beta*(1+(r/rt)**beta)**(-2.0) )+ rho_outer(r,p1) * se/r**2. *(1 + se)


def logdev2_einasto_profile(rr, p):
    p0 = p[0:3]
    p1 = p[3:5]
    p2 = p[5:8]
    rhos, alpha, rs, rho0, se, rt, beta, gamma = p
    rhos=10.0**rhos
    #rho0=10.0**rho0
    alpha = 10.0**alpha
    rs = 10.0**rs
    gamma=10.0**gamma
    beta=10.0**beta
    rt=10.0**rt
    result=np.zeros(0)
    for r in rr:
        result=np.append(result,r**2/einasto_profile(r,p)**2*(dev_einasto_profile(r,p)*einasto_profile(r,p)/r-(dev_einasto_profile(r,p))**2.0+einasto_profile(r,p)*dev2_einasto_profile(r,p)))
    return result

def dsigma(z, r, p):
    res = 2.0 * einasto_profile((z**2+r**2)**0.5, p)
    if np.isinf(res).any():
        res = 0.0
    if np.isnan(res).any():
        res = 0.0
    return res

def dev_dsigma(z, r, p):
    res = 2.0 * dev_einasto_profile((z**2+r**2)**0.5, p)*r/(z**2+r**2)**0.5
    if np.isinf(res).any(): 
        res = 0.0
    if np.isnan(res).any():
        res = 0.0
    return res

def dev2_dsigma(z, r, p):
    res = 2.0 * ( dev_einasto_profile((z**2+r**2)**0.5, p)/(z**2+r**2)**0.5 -dev_einasto_profile((z**2+r**2)**0.5, p)*r/(z**2+r**2)**(3.0/2.0)+r**2/(z**2+r**2)*dev2_einasto_profile((z**2+r**2)**0.5,p))
    if np.isinf(res).any():
        res = 0.0
    if np.isnan(res).any():
        res = 0.0
    return res

def dsigma_inner(z, r, p):
    res = 2.0 * einasto_profile_inner((z**2+r**2)**0.5, p)
    if np.isinf(res):
        print("Is inf")
        res = 0.0
    if np.isnan(res):
        print("Is nan")
        res = 0.0
    return res

def dsigma_outer(z, r, p):
    res = 2.0 * einasto_profile_outer((z**2+r**2)**0.5, p)
    if np.isinf(res):
        print("Is inf")
        res = 0.0
    if np.isnan(res):
        print("Is nan")
        res = 0.0
    return res


def rayleigh_distribution(r,sigma):
    return r/sigma**2*np.exp(-r**2/(2*sigma**2))


def surf_den(rr, p,zmax=40.0,part=0):
    '''part=0 : use full profile
       part=-1 : use inner profile
       part=1 : use outer profile
       2D density calculated from 3D einasto profile 
    '''
    result=np.zeros(0)
    if rr.size!=1:
        for r in rr:
            if part==0:
	        res = integrate.quad(dsigma, 0.0, zmax, args=(r, p))[0]
            elif part==1:
	        res = integrate.quad(dsigma_outer, 0, zmax, args=(r, p))[0]   
            else:
	        res = integrate.quad(dsigma_inner, 0, zmax, args=(r, p))[0]   
            result = np.append(result,res)
    else:
	if part==0:
	    result = integrate.quad(dsigma, 0, zmax, args=(rr, p))[0]   
	elif part==1:
	    result = integrate.quad(dsigma_outer, 0, zmax, args=(rr, p))[0]   
	else:
	    result = integrate.quad(dsigma_inner, 0, zmax, args=(rr, p))[0]   
    return result


def d_cond_surf_den(theta,rmis,r,p,part=0):
    res =  1/(2*np.pi)*surf_den(np.sqrt(r**2+rmis**2-2*r*rmis*np.cos(theta)),p[0:8],zmax=40.0,part=part)
    if np.isinf(res):
        print("Is inf")
        res = 0.0
    if np.isnan(res):
        print("Is nan")
        res = 0.0
    return res


def cond_surf_den(rmis,r,p,part=0):
    res = integrate.quad(d_cond_surf_den, 0, 2*np.pi, args=(rmis,r,p,part))[0]
    return res


def d_surf_den_mis(rmis,r,p,part=0):
    sigma=p[-1]
    res = cond_surf_den(rmis,r,p,zmax=40.0,part=part)*rayleigh_distribution(rmis,sigma)    
    if np.isinf(res):
        print("Is inf")
        res = 0.0
    if np.isnan(res):
        print("Is nan")
        res = 0.0
    return res


def d_surf_den_mis_alt(z,rmis,theta,r,p,part=0):
    sigma=p[-1]
    pmin=p[0:8]
    if part==0:
        res = 2.0*rayleigh_distribution(rmis,sigma)/(2*np.pi)*2.0*einasto_profile((z**2+r**2+rmis**2-2*r*rmis*np.cos(theta))**0.5,pmin)
    elif part==-1:
        res = 2.0*rayleigh_distribution(rmis,sigma)/(2*np.pi)*2.0*einasto_profile_inner((z**2+r**2+rmis**2-2*r*rmis*np.cos(theta))**0.5,pmin)
    else:        
        res = 2.0*rayleigh_distribution(rmis,sigma)/(2*np.pi)*2.0*einasto_profile_outer((z**2+r**2+rmis**2-2*r*rmis*np.cos(theta))**0.5,pmin)
    if np.isinf(res):
        print("Is inf")
        res = 0.0
    if np.isnan(res):
        print("Is nan")
        res = 0.0
    return res


def integrator_3D(func,x_borders=(0,0),y_borders=(0,0),z_borders=(0,0),divisions=50,args=()):
    x_range = np.linspace(x_borders[0],x_borders[1],divisions)
    x_centers = (x_range[1:]-x_range[:-1])/2.0 + x_range[:-1]
    y_range = np.linspace(y_borders[0],y_borders[1],divisions)
    y_centers = (y_range[1:]-y_range[:-1])/2.0 + y_range[:-1]
    z_range = np.linspace(z_borders[0],z_borders[1],divisions)
    z_centers = (z_range[1:]-z_range[:-1])/2.0 + z_range[:-1]
    volume = (x_range[1]-x_range[0])*(y_range[1]-y_range[0])*(z_range[1]-z_range[0])
    full = np.zeros((1,3))
    fullx = np.repeat(x_centers,divisions**2)
    foo = np.repeat(y_centers,divisions)
    fully = np.tile(foo,divisions)
    fullz = np.tile(z_centers,divisions**2)
    integral = 0.0
    for i in range(fullx.size):
       integral += func(fullx[i],fully[i],fullz[i],*args)*volume
    return integral


def spline_grid_func(R,Rmis,theta,p,zmax,part):
    sigma = p[-1]
    psub = p[0:8]
    result = rayleigh_distribution(Rmis,sigma)/(2.*np.pi)*surf_den(np.sqrt(R**2+Rmis**2+2.*R*Rmis*np.cos(theta)),psub,zmax,part)
    return result


def surf_den_mis_spline(rr,p,rmis_max=2.0,zmax=40.0,part=0,rmis_divisions=10,theta_divisions=10):
    rmis_array = np.linspace(0.0,rmis_max,rmis_divisions)
    theta_array = np.linspace(0.0,2*np.pi,theta_divisions)
    approx_array = np.zeros((rmis_divisions,theta_divisions))
    if rr.size!=1:
	result = np.zeros(0)
	for r in rr:
	    for i in range(rmis_divisions):
		for j in range(theta_divisions): 
		    approx_array[i,j] = spline_grid_func(r,rmis_array[i],theta_array[j],p,zmax,part)
	    spline = spl.RectBivariateSpline(rmis_array,theta_array,approx_array)
	    res = spline.integral(0.0,rmis_max,0.0,2.*np.pi)
	    result = np.append(result,res)
    else:
	for i in range(rmis_divisions):
	    for j in range(theta_divisions): 
		approx_array[i,j] = spline_grid_func(rr,rmis_array[i],theta_array[j],p,zmax,part)
	spline = spl.RectBivariateSpline(rmis_array,theta_array,approx_array)
	resuls = spline.integral(0.0,rmis_max,0.0,2.*np.pi)
    return result


def surf_den_mis_spline_alt(rr,p,rmis_max=2.0,zmax=40.0,part=0,rmis_divisions=10,theta_divisions=10):
    rmis_array = np.linspace(0.0,rmis_max,rmis_divisions)[:-1]
    theta_array = np.linspace(0.0,2*np.pi,theta_divisions)
    binsize=rmis_max*2.*np.pi/(theta_divisions-1.)
    result = np.zeros(0)
    for r in rr:
	approx_array_2=np.zeros(rmis_divisions-1)
	for i in range(rmis_divisions-1):
	    approx_array_1=np.zeros(theta_divisions)
	    for j in range(theta_divisions): 
		approx_array_1[j] = spline_grid_func(r,rmis_array[i],theta_array[j],p,zmax,part)
	    spline_1 = spl.UnivariateSpline(theta_array,approx_array_1)
	    approx_array_2[i]=spline_1.integral(0.0,2.*np.pi)
	res = np.sum(approx_array_2)*binsize
	result = np.append(result,res)
    return result


def surf_den_mis_alt(rr,p,rmis_max=2.0,zmax=40.0,part=0):
    result = np.zeros(0)
    if rr.size!=1:
        for r in rr:
	    res = integrate.tplquad(d_surf_den_mis_alt,0.,np.pi,lambda theta:0.0,lambda theta:rmis_max,lambda theta,rmis:0.0,lambda theta,rmis:zmax,args=(r,p,part),epsrel=1e-2)[0]
	    result = np.append(result,res)
    else:
	result = integrate.tplquad(d_surf_den_mis_alt,0.,np.pi,lambda theta:0.0,lambda theta:rmis_max,lambda theta,rmis:0.0,lambda theta,rmis:zmax,args=(rr,p,part),epsrel=1e-2)[0]
    return result


def surf_den_mis_primitive(rr,p,rmis_max=2.0,zmax=40.0,part=0):
    result = np.zeros(0)
    if rr.size!=1:
        for r in rr:
	    res = integrator_3D(d_surf_den_mis_alt,(0.0,zmax),(0.0,rmis_max),(0.0,np.pi),divisions=100,args=(r,p,part))
	    result = np.append(result,res)
    else:
	    result = integrator_3D(d_surf_den_mis_alt,(0.0,zmax),(0.0,rmis_max),(0.0,np.pi),divisions=100,args=(rr,p,part))
    return result


def surf_den_mis(rr,p,rmis_max=2.0,part=0):
    result = np.zeros(0)
    if rr.size!=1:
        for r in rr:
            res = integrate.quad(d_surf_den_mis, 0, rmis_max, args=(r, p, part))[0]
            result = np.append(result,res)
    else:
        result = integrate.quad(d_surf_den_mis, 0, rmis_max, args=(rr, p, part))[0]
    return result


def total_surf_den(rr,p,zmax=40.0,part=0,rmis_divisions=5,theta_divisions=5,alt=0):
    if alt==0:
        result = (1. - p[8])*surf_den(rr,p[0:8],part=part) + p[8]*surf_den_mis_spline(rr,p,part=part,rmis_divisions=rmis_divisions,theta_divisions=theta_divisions)
    else:
        result = (1. - p[8])*surf_den(rr,p[0:8],part=part) + p[8]*surf_den_mis_spline_alt(rr,p,part=part,rmis_divisions=rmis_divisions,theta_divisions=theta_divisions)
    return result


def dev_surf_den(rr, p,zmax=40.0):
    '''derivative of 2D density calculated from 3D einasto profile '''
    result=np.zeros(0)
    if rr.size!=1:
        for r in rr:
	    res = integrate.quad(dev_dsigma, 0, zmax, args=(r, p))[0]   
            result = np.append(result,res)
    else:
	result=integrate.quad(dev_dsigma, 0, zmax, args=(rr, p))[0]
    return result


def logdev_surf_den(rr, p, recursive=False, surf=[],zmax=40.0):
    '''derivative of 2D density calculated from 3D einasto profile '''
    result=np.zeros(0)
    if rr.size==1:
	return rr / integrate.quad(dsigma, 0, zmax, args=(rr, p))[0] * integrate.quad(dev_dsigma, 0, zmax, args=(rr, p))[0]   
    for i,r in enumerate(rr):
        if recursive==False:
	    res = r / integrate.quad(dsigma, 0, zmax, args=(r, p))[0] * integrate.quad(dev_dsigma, 0, zmax, args=(r, p))[0]   
        else:
	    res =  r / surf[i] * integrate.quad(dev_dsigma, 0, zmax, args=(r, p))[0]   
        result = np.append(result,res)
    return result


def dev2_surf_den(rr, p,zmax=40.0):
    '''second derivative of 2D density calculated from 3D einasto profile '''
    result=np.zeros(0)
    if rr.size!=1:
        for r in rr:
	    res = integrate.quad(dev2_dsigma, 0, zmax, args=(r, p))[0]   
            result = np.append(result,res)
    else:
	result=integrate.quad(dev2_dsigma, 0, zmax, args=(rr, p))[0]
    return result

def logdev2_surf_den(rr, p, recursive=False,surf=[],surf_dev=[],surf_dev2=[],zmax=40.0):
    '''second derivative of 2D density calculated from 3D einasto profile '''
    result=np.zeros(0)
    for i,r in enumerate(rr):
        if recursive==False:
            res =r**2/(integrate.quad(dsigma, 0, zmax, args=(r, p))[0])**2*(integrate.quad(dsigma, 0, zmax, args=(r, p))[0]*integrate.quad(dev_dsigma, 0, zmax, args=(r, p))[0]/r-(integrate.quad(dev_dsigma, 0, zmax, args=(r, p))[0])**2+integrate.quad(dsigma, 0, zmax, args=(r, p))[0]*integrate.quad(dev2_dsigma, 0, zmax, args=(r, p))[0])   
        else:
            res=r**2/surf[i]**2*(surf[i]*surf_dev[i]/r-surf_dev[i]**2+surf[i]*surf_dev2[i])
        result = np.append(result,res)
    return result


def surf_den_inner(rr, p,zmax=40.0):
    '''2D density calculated from 3D einasto profile '''
    result=np.zeros(0)
    for r in rr:
	res = integrate.quad(dsigma_inner, 0, zmax, args=(r, p))[0]
        result = np.append(result,res)
    return result

def surf_den_outer(rr, p, mc, zmax=40.0):
    '''2D density calculated from 3D einasto profile '''
    result=np.zeros(0)
    for r in rr:
	res = integrate.quad(dsigma_outer, 0, zmax, args=(r, p,mc))[0]
        result = np.append(result,res)
    return result

def get_chisquare_reload(p, rr, yy, yy_err, invcov, prior,type_="",mc=False,alt=False):
    """Calculates sum of chisquares between a model function given by the parameters p and some data points given by yy, yy_err and the covariance matrix of those (cov). """
    if mc==False:
        rhos, alpha, rs, rho0, se, rt, beta, gamma =p
        yymodel = surf_den(rr, p)
    else:
        rhos, alpha, rs, rho0, se, rt, beta, gamma, fmin, sigma =p
        yymodel = total_surf_den(rr, p,alt=alt)
    if (type_=="Planck_PS_21_deproj") | (type_=="Planck_PS_21.5_deproj") | (type_=="Planck_PS_22_deproj"):
        yymodel=einasto_profile(rr,p)
    chi_sq=np.transpose(yymodel-yy).dot(invcov).dot(yymodel-yy)
    chisq_prior=0

    if prior=="best":
        #Double prior size for alpha,beta,gamma
        chisq_prior += ((beta-np.log10(6.0))/0.4)**2
        chisq_prior += ((gamma-np.log10(4.0))/0.4)**2
        chisq_prior += ((alpha-np.log10(0.2))/1.2)**2
        if mc==True:
	    chisq_prior += ((fmin-0.15)/0.21)**2
	    chisq_prior += ((sigma-0.41)/0.3)**2
        if 10**rs>5.0:
            chisq_prior += ((10**rs-5.0)/0.00001)**2
        if 10**rs<0.1:
            chisq_prior += ((10**rs-0.1)/0.00001)**2
        if 10**rt>5.0:
            chisq_prior += ((10**rt-5.0)/0.00001)**2
        if 10**rt<0.1:
            chisq_prior += ((10**rt-0.1)/0.00001)**2
        if se<0.0:
            chisq_prior+=((se-0.0)/0.00001)**2
        if rho0<0.0:
            chisq_prior+=((rho0-0.0)/0.00001)**2
        #if se>1.8:
        #    chisq_prior+=((se-1.8)/0.00001)**2
        
    if prior=="None":
        pass


    chi_sq+=chisq_prior
    return chi_sq
def get_chisquare(p, rr, yy, yy_err, invcov, prior,type_="",est="",mc=False):
    """Calculates sum of chisquares between a model function given by the parameters p and some data points given by yy, yy_err and the covariance matrix of those (cov). """
    if mc==False:
        rhos, alpha, rs, rho0, se, rt, beta, gamma =p
    else:
        rhos, alpha, rs, rho0, se, rt, beta, gamma, fmin, sigma =p

    yymodel = surf_den(rr, p,mc)
    if (type_=="Planck_PS_21_deproj") | (type_=="Planck_PS_21.5_deproj") | (type_=="Planck_PS_22_deproj"):
        yymodel=einasto_profile(rr,p)
    chi_sq=np.transpose(yymodel-yy).dot(invcov).dot(yymodel-yy)
    chisq_prior=0

    if prior=="surhud":
        chisq_prior += ((beta-np.log10(6.0))/0.2)**2
        chisq_prior += ((gamma-np.log10(4.0))/0.2)**2
        chisq_prior += ((alpha-np.log10(0.2))/0.6)**2
        if rs>5.0:
            chisq_prior += ((10**rs-5.0)/0.00001)**2
        if rs<0.1:
            chisq_prior += ((10**rs-0.1)/0.00001)**2
        if rt>5.0:
            chisq_prior += ((10**rt-5.0)/0.00001)**2
        if rt<0.1:
            chisq_prior += ((10**rt-0.1)/0.00001)**2
    if prior=="self":
        chisq_prior += ((beta-0.8)/0.6)**2
        chisq_prior += ((gamma-0.5)/0.6)**2
        chisq_prior += ((alpha-(-np.log10(0.2)))/0.6)**2
        if rs>5.0:
            chisq_prior += ((rs-5.0)/0.00001)**2
        if rs<0.1:
            chisq_prior += ((rs-0.1)/0.00001)**2
        if rt>5.0:
            chisq_prior += ((rt-5.0)/0.00001)**2
        if rt<0.1:
             chisq_prior += ((rt-0.1)/0.00001)**2
    if prior=="MLE":
        mlep = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro/Mest/mles/MLE_parameters_"+str(type_)+"_("+str(est)+").txt", header=None)
        mlep = np.asarray(mlep.values).reshape(mlep.size)
        rhos_mle=mlep[0]
        alpha_mle=mlep[1]
        rs_mle=mlep[2]
        rho0_mle=mlep[3]
        se_mle=mlep[4]
        rt_mle=mlep[5]
        beta_mle=mlep[6]
        gamma_mle=mlep[7]
        chisq_prior += ((beta-beta_mle)/0.6)**2
        chisq_prior += ((gamma-gamma_mle)/0.6)**2
        chisq_prior += ((alpha-alpha_mle)/0.6)**2
         
        rs_max=rs_mle+0.25
        rs_min=rs_mle-0.25
        rt_min=rt_mle-0.25
        rt_max=rt_mle+0.25
        if rs>rs_max:
            chisq_prior += ((rs-rs_max)/0.00001)**2
        if rs<rs_min:
            chisq_prior += ((rs-rs_min)/0.00001)**2
        if rt>rs_max:
            chisq_prior += ((rt-rt_max)/0.00001)**2
        if rt<rs_min:
             chisq_prior += ((rt-rt_min)/0.00001)**2

    if prior=="long":
        mlep = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro/Mest/mles/MLE_parameters_"+str(type_)+"_("+str(est)+").txt", header=None)
        mlep = np.asarray(mlep.values).reshape(mlep.size)
        rhos_mle=mlep[0]
        alpha_mle=mlep[1]
        rs_mle=mlep[2]
        rho0_mle=mlep[3]
        se_mle=mlep[4]
        rt_mle=mlep[5]
        beta_mle=mlep[6]
        gamma_mle=mlep[7]

        rhos_sigma=2.2
        alpha_sigma=0.9
        rs_sigma=1.3
        rho0_sigma=0.024
        se_sigma=0.35
        rt_sigma=0.18
        beta_sigma=0.61
        gamma_sigma=0.61

        over=1.2

        chisq_prior += ((rhos-beta_mle)/(rhos_sigma*over))**2
        chisq_prior += ((alpha-beta_mle)/(alpha_sigma*over))**2
        chisq_prior += ((rs-beta_mle)/(rs_sigma*over))**2
        chisq_prior += ((rho0-beta_mle)/(rho0_sigma*over))**2
        chisq_prior += ((se-beta_mle)/(se_sigma*over))**2
        chisq_prior += ((rt-beta_mle)/(rt_sigma*over))**2
        chisq_prior += ((beta-beta_mle)/(beta_sigma*over))**2
        chisq_prior += ((gamma-beta_mle)/(gamma_sigma*over))**2
         

    if prior=="flat":
        mlep = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro/Mest/mles/MLE_parameters_"+str(type_)+"_("+str(est)+").txt", header=None)
        mlep = np.asarray(mlep.values).reshape(mlep.size)
        rhos_mle=mlep[0]
        alpha_mle=mlep[1]
        rs_mle=mlep[2]
        rho0_mle=mlep[3]
        se_mle=mlep[4]
        rt_mle=mlep[5]
        beta_mle=mlep[6]
        gamma_mle=mlep[7]

        snr_red=263.0
        snr_sz=60.0
        incr=snr_red/snr_sz

        rhos_sigma=0.4/0.6*incr*abs(rhos_mle)
        alpha_sigma=0.2/0.98*incr*abs(alpha_mle)
        rs_sigma=0.27/0.32*incr*abs(rs_mle)
        rho0_sigma=0.055/0.055*incr*abs(rho0_mle)
        se_sigma=0.08/1.55*incr*abs(se_mle)
        rt_sigma=0.04/0.087*incr*abs(rt_mle)
        beta_sigma=0.1/1.0*incr*abs(beta_mle)
        gamma_sigma=0.14/0.85*incr*abs(gamma_mle)

        over=3.0 #disregard everything that is more than 3 sigma away from MLE result

        rhos_min=rhos_mle-over*rhos_sigma
        rhos_max=rhos_mle+over*rhos_sigma
        alpha_min=alpha_mle-over*alpha_sigma
        alpha_max=alpha_mle+over*alpha_sigma
        rs_min=rs_mle-over*rs_sigma
        rs_max=rs_mle+over*rs_sigma
        rho0_min=rho0_mle-over*rho0_sigma
        rho0_max=rho0_mle+over*rho0_sigma

        #Putting stronger constrain on the problematic se parameter
        se_min=se_mle-over*se_sigma
        se_max=se_mle+over*se_sigma

        rt_min=rt_mle-over*rt_sigma
        rt_max=rt_mle+over*rt_sigma
        beta_min=beta_mle-over*beta_sigma
        beta_max=beta_mle+over*beta_sigma
        gamma_min=gamma_mle-over*gamma_sigma
        gamma_max=gamma_mle+over*gamma_sigma
 
        #print("Alpha flat prior from: "+str(alpha_min)+" to "+str(alpha_max))
        #print("rs flat prior from: "+str(rs_min)+" to "+str(rs_max))
        #print("rt flat prior from: "+str(rt_min)+" to "+str(rt_max))
        #print("Beta flat prior from: "+str(beta_min)+" to "+str(beta_max))
        #print("Gamma flat prior from: "+str(gamma_min)+" to "+str(gamma_max))


        if rhos>rhos_max:
            chisq_prior += ((rhos-rhos_max)/0.00001)**2
        if rhos<rhos_min:
            chisq_prior += ((rhos-rhos_min)/0.00001)**2
        if alpha>alpha_max:
            chisq_prior += ((alpha-alpha_max)/0.00001)**2
        if alpha<alpha_min:
            chisq_prior += ((alpha-alpha_min)/0.00001)**2
        if rs>rs_max:
            chisq_prior += ((rs-rs_max)/0.00001)**2
        if rs<rs_min:
            chisq_prior += ((rs-rs_min)/0.00001)**2
        if rho0>rho0_max:
            chisq_prior += ((rho0-rho0_max)/0.00001)**2
        if rho0<rho0_min:
            chisq_prior += ((rho0-rho0_min)/0.00001)**2
        if se>se_max:
            chisq_prior += ((se-se_max)/0.00001)**2
        if se<se_min:
            chisq_prior += ((se-se_min)/0.00001)**2
        if rt>rt_max:
            chisq_prior += ((rt-rt_max)/0.00001)**2
        if rt<rt_min:
            chisq_prior += ((rt-rt_min)/0.00001)**2
        if beta>beta_max:
            chisq_prior += ((beta-beta_max)/0.00001)**2
        if beta<beta_min:
            chisq_prior += ((beta-beta_min)/0.00001)**2
        if gamma>gamma_max:
            chisq_prior += ((gamma-gamma_max)/0.00001)**2
        if gamma<gamma_min:
            chisq_prior += ((gamma-gamma_min)/0.00001)**2


    if prior=="flat2":
        #Reading the MLE values
        mlep = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro/Mest/mles/MLE_parameters_"+str(type_)+"_"+str(est)+".txt", header=None)
        mlep = np.asarray(mlep.values).reshape(mlep.size)
        alpha_mle=mlep[1]
        rs_mle=mlep[2]
        rt_mle=mlep[5]
        beta_mle=mlep[6]
        gamma_mle=mlep[7]

        snr_red=263.0 #SNR of RedMaPPer
        snr_sz=60.0 #SNR of Planck SZ
        incr=snr_red/snr_sz

        #get ranges for the new priors from the RedMaPPer errors
        #The structure is new_error= RedMaPPer_error * SNR_correction

        alpha_sigma=0.2*incr  # SZ_error= RedMaPPer_error/
        rs_sigma=0.27*incr
        rt_sigma=0.04*incr
        beta_sigma=0.1*incr
        gamma_sigma=0.14*incr

        over=3.0 #disregard everything that is more than 3 sigma away from MLE result

        alpha_min=alpha_mle-over*alpha_sigma
        alpha_max=alpha_mle+over*alpha_sigma
        rs_min=rs_mle-over*rs_sigma
        rs_max=rs_mle+over*rs_sigma

        over_rt=over*6.0
        rt_min=rt_mle-over_rt*rt_sigma
        rt_max=rt_mle+over_rt*rt_sigma

        beta_min=beta_mle-over*beta_sigma
        beta_max=beta_mle+over*beta_sigma
        gamma_min=gamma_mle-over*gamma_sigma
        gamma_max=gamma_mle+over*gamma_sigma

        print("Alpha flat prior from: "+str(alpha_min)+" to "+str(alpha_max))
        print("rs flat prior from: "+str(rs_min)+" to "+str(rs_max))
        print("rt flat prior from: "+str(rt_min)+" to "+str(rt_max))
        print("Beta flat prior from: "+str(beta_min)+" to "+str(beta_max))
        print("Gamma flat prior from: "+str(gamma_min)+" to "+str(gamma_max))

 
        if alpha>alpha_max:
            chisq_prior += ((alpha-alpha_max)/0.00001)**2
        if alpha<alpha_min:
            chisq_prior += ((alpha-alpha_min)/0.00001)**2
        if rs>rs_max:
            chisq_prior += ((rs-rs_max)/0.00001)**2
        if rs<rs_min:
            chisq_prior += ((rs-rs_min)/0.00001)**2
        if rt>rt_max:
            chisq_prior += ((rt-rt_max)/0.00001)**2
        if rt<rt_min:
            chisq_prior += ((rt-rt_min)/0.00001)**2
        if beta>beta_max:
            chisq_prior += ((beta-beta_max)/0.00001)**2
        if beta<beta_min:
            chisq_prior += ((beta-beta_min)/0.00001)**2
        if gamma>gamma_max:
            chisq_prior += ((gamma-gamma_max)/0.00001)**2
        if gamma<gamma_min:
            chisq_prior += ((gamma-gamma_min)/0.00001)**2

    if prior=="flat2_tight":
        #Reading the MLE values
        mlep = pd.read_csv("/home/dominik.zuercher/Documents/RSP_Pro/Mest/mles/MLE_parameters_"+str(type_)+"_"+str(est)+".txt", header=None)
        mlep = np.asarray(mlep.values).reshape(mlep.size)
        alpha_mle=mlep[1]
        rs_mle=mlep[2]
        rt_mle=mlep[5]
        beta_mle=mlep[6]
        gamma_mle=mlep[7]

        snr_red=263.0 #SNR of RedMaPPer
        snr_sz=60.0 #SNR of Planck SZ
        incr=snr_red/snr_sz

        #get ranges for the new priors from the RedMaPPer errors
        #The structure is new_error= RedMaPPer_error * SNR_correction

        alpha_sigma=0.2*incr  # SZ_error= RedMaPPer_error/
        rs_sigma=0.27*incr
        rt_sigma=0.04*incr
        beta_sigma=0.1*incr
        gamma_sigma=0.14*incr

        over=1.0 #disregard everything that is more than 3 sigma away from MLE result

        alpha_min=alpha_mle-over*alpha_sigma
        alpha_max=alpha_mle+over*alpha_sigma
        rs_min=rs_mle-over*rs_sigma
        rs_max=rs_mle+over*rs_sigma

        over_rt=over*1.0
        rt_min=rt_mle-over_rt*rt_sigma
        rt_max=rt_mle+over_rt*rt_sigma

        beta_min=beta_mle-over*beta_sigma
        beta_max=beta_mle+over*beta_sigma
        gamma_min=gamma_mle-over*gamma_sigma
        gamma_max=gamma_mle+over*gamma_sigma

        #print("Alpha flat prior from: "+str(alpha_min)+" to "+str(alpha_max))
        #print("rs flat prior from: "+str(rs_min)+" to "+str(rs_max))
        #print("rt flat prior from: "+str(rt_min)+" to "+str(rt_max))
        #print("Beta flat prior from: "+str(beta_min)+" to "+str(beta_max))
        #print("Gamma flat prior from: "+str(gamma_min)+" to "+str(gamma_max))

 
        if alpha>alpha_max:
            chisq_prior += ((alpha-alpha_max)/0.00001)**2
        if alpha<alpha_min:
            chisq_prior += ((alpha-alpha_min)/0.00001)**2
        if rs>rs_max:
            chisq_prior += ((rs-rs_max)/0.00001)**2
        if rs<rs_min:
            chisq_prior += ((rs-rs_min)/0.00001)**2
        if rt>rt_max:
            chisq_prior += ((rt-rt_max)/0.00001)**2
        if rt<rt_min:
            chisq_prior += ((rt-rt_min)/0.00001)**2
        if beta>beta_max:
            chisq_prior += ((beta-beta_max)/0.00001)**2
        if beta<beta_min:
            chisq_prior += ((beta-beta_min)/0.00001)**2
        if gamma>gamma_max:
            chisq_prior += ((gamma-gamma_max)/0.00001)**2
        if gamma<gamma_min:
            chisq_prior += ((gamma-gamma_min)/0.00001)**2

    if prior=="best":
        #Double prior size for alpha,beta,gamma
        chisq_prior += ((beta-np.log10(6.0))/0.4)**2
        chisq_prior += ((gamma-np.log10(4.0))/0.4)**2
        chisq_prior += ((alpha-np.log10(0.2))/1.2)**2
        if mc==True:
	    chisq_prior += ((fmin-0.2)/0.07)**2
	    chisq_prior += ((sigma-0.392)/0.1)**2
        if 10**rs>5.0:
            chisq_prior += ((10**rs-5.0)/0.00001)**2
        if 10**rs<0.1:
            chisq_prior += ((10**rs-0.1)/0.00001)**2
        if 10**rt>5.0:
            chisq_prior += ((10**rt-5.0)/0.00001)**2
        if 10**rt<0.1:
            chisq_prior += ((10**rt-0.1)/0.00001)**2
        if se<0.0:
            chisq_prior+=((se-0.0)/0.00001)**2
        if rho0<0.0:
            chisq_prior+=((rho0-0.0)/0.00001)**2
        #if se>1.8:
        #    chisq_prior+=((se-1.8)/0.00001)**2
        
    if prior=="None":
        pass


    chi_sq+=chisq_prior
    return chi_sq


def chisquare_opt(p, *args):
    """Chisquare function suitable for optimizing """
    chisq=get_chisquare(p, *args)
    if np.isnan(chisq):
        print("Nan")
        return -np.inf
    if not np.isfinite(chisq):
        print("inf")
        return -np.inf
    return 0.5*chisq

def logdev2d_opt(r, *args):
    """Chisquare function suitable for optimizing """
    logdev=logdev_surf_den(r,p=args)
    return logdev

def mindev2_procedure(p):
    rmin = fminbound(logdev2d_opt, 0.0,5.0, args=(p))
    return rmin

def logdev3d_opt(r, *args):
    """Chisquare function suitable for optimizing """
    logdev=logdev_einasto_profile(r, p=args)
    return logdev

def mindev3_procedure(p):
    rmin = fminbound(logdev3d_opt,0.0,9.0, args=(p))
    return rmin


def MLE_procedure(rr, yy, yy_err, invcov, p0, prior,type_,est):
    """p0 is inital paramater set (8 values)"""
    pmle = fmin(chisquare_opt, p0, args=(rr, yy, yy_err, invcov, prior,type_,est))
    return pmle

def lnlike(p, rr, yy, yy_err, invcov, prior,type_,est):
    chisq = get_chisquare(p, rr, yy, yy_err, invcov, prior,type_=type_,est=est)
    if np.isnan(chisq):
        print("Nan")
        return -np.inf
    if not np.isfinite(chisq):
        print("inf")
        return -np.inf
    return -0.5*chisq

def lnlike_reload(p, rr, yy, yy_err, invcov, prior,type_,mc,alt=False):
    chisq = get_chisquare_reload(p, rr, yy, yy_err, invcov, prior,type_,mc,alt)
    if np.isnan(chisq):
        print("Nan")
        return -np.inf
    if not np.isfinite(chisq):
        print("inf")
        return -np.inf
    return -0.5*chisq

def MLE_procedure_alt(rr, yy, yy_err, invcov, p0, prior,type_,est):
    """p0 is inital paramater set (8 values)"""
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, p0, args=(rr, yy, yy_err,invcov,prior,type_,est))
    return result["x"]

def MLE_procedure_alt_reload(rr, yy, yy_err, invcov, p0, prior,type_,mc,alt=False):
    """p0 is inital paramater set (8 values)"""
    nll = lambda *args: -lnlike_reload(*args)
    result = op.minimize(nll, p0, args=(rr, yy, yy_err,invcov,prior,type_,mc,alt))
    return result["x"]

def deproj_lnlike(p, rr, yy, yy_err, invcov, prior,type_):
    chisq = get_chisquare(p, rr, yy, yy_err, invcov, prior,type_=type_)
    if np.isnan(chisq):
        print("Nan")
        return -np.inf
    if not np.isfinite(chisq):
        print("inf")
        return -np.inf
    return -0.5*chisq

def deproj_MLE_procedure_alt(rr, yy, yy_err, invcov, p0, prior,type_):
    """p0 is inital paramater set (8 values)"""
    nll = lambda *args: -deproj_lnlike(*args)
    result = op.minimize(nll, p0, args=(rr, yy, yy_err,invcov,prior,type_))
    return result["x"]

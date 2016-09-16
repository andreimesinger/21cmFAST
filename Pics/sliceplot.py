#!/usr/bin/env python
# SIMPLEST USAGE: python sliceplot.py -i file1 file2 file3...
#
# More complex options:
USAGE = "USAGE: python sliceplot.py [--savearray] [--zindex=<z-index of slice>] [--delzindex=<offset for subtracting image>] [--filter=<smoothing sigma>] [--filterx=<x-axis smoothing sigma>] [--filtery=<y-axis smoothing sigma>] [--filterz=<z-axis smoothing sigma>] [--min=<min of plot>] [--max=<max of plot>] -i <filename1> <filename2>..."

###### LIST OF OPTIONAL ARGUMENTS:
###  --savearray is a flag indicating you want to also save the 2D slice as a .npy file.
###  --zindex= lets you specify which array index cut through the z axis (usualy the LOS axis).  DEFAULT is the midpoint, i.e. DIM/2
###  --delzindex= if this is specified, then we will plot the difference between slices at array[:,:,zindex] - array[:,:,zindex+delzindex]
###  --filter= allows you to smooth the array with a Gaussian filter with the specified standard deviation (in units of array cells).  DEFAULT is no smoothing.
###  --filterx= smooth only along the horizontal axis.  DEFAULT is no smoothing.
###  --filtery= smooth only along the vertical axis.  DEFAULT is no smoothing.
###  --filterz= smooth only along the line of sight axis.  DEFAULT is no smoothing.
###  --min= specify a minimum value for the plotted range
###  --max= specify a maximum value for the plotted range

from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import scipy
from scipy import ndimage
from matplotlib.ticker import *
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from os.path import basename
import os
import sys, getopt


#To normalize the midpoint of the colorbar
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def load_binary_data(filename, dtype=np.float32): 
     """ 
     We assume that the data was written 
     with write_binary_data() (little endian). 
     """ 
     f = open(filename, "rb") 
     data = f.read() 
     f.close() 
     _data = np.fromstring(data, dtype) 
     if sys.byteorder == 'big':
       _data = _data.byteswap()
     return _data 


#default arguements
x_sigma = -1 # if negative, then do not smooth the field
y_sigma = -1
z_sigma = -1
iso_sigma = -1
files_in = []
z_index = -1
minrange = 1e5
maxrange = -1e5
savefile = 0
del_z_index = int(0)

# check for optional arguments
if(1):
  try:
    opts, args = getopt.getopt(sys.argv[1:],"u:f:x:z:y:i:", ["filter=", "filterx=", "filtery=", "filterz=", "fx=", "fy=", "fz=", "zindex=", "min=", "max=", "savearray", "delzindex="])
  except getopt.GetoptError:
    print USAGE
    sys.exit(2)
  for opt, arg in opts:
#    print opt,arg
    if opt in ("-u", "--u", "-h", "--h", "--help"):
      print USAGE
      sys.exit()
    elif opt in ("-x", "-filterx", "--filterx"):
      x_sigma = float(arg)
    elif opt in ("-y", "-filtery", "--filtery"):
      y_sigma = float(arg)
    elif opt in ("-z", "-filterz", "--filterz"):
      z_sigma = float(arg)
    elif opt in ("-f", "--f", "-filter", "--filter"):
      iso_sigma = float(arg)
    elif opt in ("-i", "--i"):
      files_in = arg
      files_in = files_in.split()
    elif opt in ("-zindex", "--zindex"):
      z_index = int(arg)
    elif opt in ("-delzindex", "--delzindex"):
      del_z_index = int(arg)
    elif opt in ("-min", "--min"):
      minrange = float(arg)
    elif opt in ("-max", "--max"):
      maxrange = float(arg)
    elif opt in ("--savearray"):
      savefile = 1

if not files_in:
    print "No files for processing... Have you included a '-i' flag before the filenames?\n"+USAGE
    sys.exit()
      
# go through list of files and process each one
for path in files_in:
    #path = sys.argv[i]
    print 'Processing input file:'
    print '  '+path
    filename="" + path.split("/")[-1]
    
    # lightcone?
    if basename(filename)[-11:]=='lighttravel':
        DIM = int("" + path.split("_")[-3])
        label=str("" + path.split("_")[-2])

    else:
        DIM = int("" + path.split("_")[-2])
        label=str("" + path.split("_")[-1])

    if z_index < 0:
        z_index = DIM/2

    # read in the data cube located in 21cmFast/Boxes/delta_T*
    data1 = load_binary_data(path)
    data1.shape = (DIM, DIM, DIM)
    data1 = data1.reshape((DIM, DIM, DIM), order='F')

    # smooth the field?
    if iso_sigma > 0:
        print "Smoothing the entire cube with a Gassian filter of width="+str(iso_sigma)
        data1 = scipy.ndimage.filters.gaussian_filter(data1, sigma=iso_sigma)
    else:
        if x_sigma > 0:
            print "Smoothing along the x (horizontal) axis with a Gassian filter of width="+str(x_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=x_sigma, axis=1)
        if y_sigma > 0:
            print "Smoothing along the y (vertical) axis with a Gassian filter of width="+str(y_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=y_sigma, axis=0)
        if z_sigma > 0:
            print "Smoothing along the z (line of sight) axis with a Gassian filter of width="+str(z_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=z_sigma, axis=2)

        
    fig = plt.figure(dpi=72)
    sub_fig = fig.add_subplot(111)
    print "Taking a slice along the LOS direction at index="+str(z_index)
    slice = data1[:,:,z_index]

    if del_z_index: #difference image is wanted
        other_z_index = int(z_index+del_z_index)
        print "Subtracting the slice at index="+str(other_z_index)
        slice = slice - data1[:,:,other_z_index]

    
    # check box type to determine default plotting options
    # check if it is a 21cm brightness temperature box
    if basename(filename)[0:3]=='del':
        if minrange > 1e4:
            minrange = -210
        if maxrange < -1e4:
            maxrange = 30
        cmap = LinearSegmentedColormap.from_list('mycmap', ['yellow','red','black','green','blue'])
        norm = MidpointNormalize(midpoint=0)
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
        c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm \delta T_b [\mathrm{mK}]}$', fontsize=24, rotation=-90, labelpad=32)
        tick_array = np.linspace(minrange, maxrange, 8)

    # check if it is a neutral fraction box
    elif basename(filename)[0:3]=='xH_':
        if minrange > 1e4:
            minrange = 0
        if maxrange < -1e4:
            maxrange = 1
        cmap = LinearSegmentedColormap.from_list('mycmap', ['white','black'])
        norm = MidpointNormalize(midpoint=0.5)
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
        c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm x_{HI}}$', fontsize=24, rotation=-90, labelpad=32)
        tick_array = np.linspace(minrange, maxrange, 6)

    # check it is a density box
    elif basename(filename)[0:3]=='upd':
        if minrange > 1e4:
            minrange = -0.5
        if maxrange < -1e4:
            maxrange = 0.5
        slice = np.log10(1+data1[:,:,z_index])
        cmap = LinearSegmentedColormap.from_list('mycmap', ['darkblue', 'green', 'yellow', 'red'])
        norm = MidpointNormalize(midpoint=0)
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
        c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm log(\Delta)}$', fontsize=24, rotation=-90, labelpad=32)
        tick_array = np.linspace(minrange, maxrange, 5)


    c_bar.set_ticks(tick_array)
	
    for t in c_bar.ax.get_yticklabels():
        t.set_fontsize(14)

    if del_z_index:
        endstr = '_zindex'+str(z_index)+'-'+str(z_index+del_z_index)
    else:
        endstr = '_zindex'+str(z_index)
    plt.savefig(filename+endstr+'.pdf', bbox_inches='tight')

    # do we want to save the array file?
    if savefile:
        np.save(filename+endstr, slice)

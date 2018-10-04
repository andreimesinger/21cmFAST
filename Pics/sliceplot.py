#!/usr/bin/env python
# SIMPLEST USAGE: python sliceplot.py -i file1 file2 file3...
#
#  Default slice is along the x-direction at DIM/2
#
# More complex options:
USAGE = "USAGE: python sliceplot.py [--savearray] [--zindex=<z-index of slice>] [--delzindex=<offset for subtracting image>] [--filter=<smoothing sigma>] [--filterx=<x-axis smoothing sigma>] [--filtery=<y-axis smoothing sigma>] [--filterz=<z-axis smoothing sigma>] [--min=<min of plot>] [--max=<max of plot>] -i <filename1> <filename2>..."

###### LIST OF OPTIONAL ARGUMENTS:
###  --savearray is a flag indicating you want to also save the 2D slice as a .npy file.
###  --zindex= lets you specify which array index cut through the z axis (usualy the LOS axis).
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
import sys, argparse


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

# Parse the command line options
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input filenames', nargs='+', required=True)
parser.add_argument("-f", "--filter", type=float, default=-1, help="smooth the array with a Gaussian filter with the specified standard deviation (in units of array cells).  DEFAULT is no smoothing.")
parser.add_argument("-x", "--filterx", type=float, default=-1, help="smooth only along the horizontal axis.  DEFAULT is no smoothing.")
parser.add_argument("-y", "--filtery", type=float, default=-1, help="smooth only along the vertical axis.  DEFAULT is no smoothing.")
parser.add_argument("-z", "--filterz", type=float, default=-1, help="smooth only along the LOS (z) axis.  DEFAULT is no smoothing.")
parser.add_argument("--zindex", type=int, default=-1, help="specify which array index cut through the z axis (usualy the LOS axis).  DEFAULT is the midpoint, i.e. DIM/2.")
parser.add_argument("--delzindex", type=int, default=-1, help="if this is specified, then we will plot the difference between slices at array[:,:,zindex] - array[:,:,zindex+delzindex]")
parser.add_argument("--min", type=float, default=1e5, help="specify a minimum value for the plotting range.")
parser.add_argument("--max", type=float, default=-1e5, help="specify a maximum value for the plotting range.")
parser.add_argument("--savearray", help="flag indicating you want to also save the 2D slice as a .npy file.", action="store_true")

args = parser.parse_args()

# go through list of files and process each one
for path in args.input:

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

    if args.zindex >= 0:
        z_index = args.zindex

    # read in the data cube located in 21cmFast/Boxes/delta_T*
    data1 = load_binary_data(path)
    data1.shape = (DIM, DIM, DIM)
    data1 = data1.reshape((DIM, DIM, DIM), order='F')

    # smooth the field?
    if args.filter >= 0:
        iso_sigma = args.filter
        print "Smoothing the entire cube with a Gassian filter of width="+str(iso_sigma)
        data1 = scipy.ndimage.filters.gaussian_filter(data1, sigma=iso_sigma)
    else:
        if args.filterx >= 0:
            x_sigma = args.filterx
            print "Smoothing along the x (horizontal) axis with a Gassian filter of width="+str(x_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=x_sigma, axis=1)
        if args.filtery >= 0:
            y_sigma = args.filtery
            print "Smoothing along the y (vertical) axis with a Gassian filter of width="+str(y_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=y_sigma, axis=0)
        if args.filterz >= 0:
            z_sigma = args.filterz
            print "Smoothing along the z (line of sight) axis with a Gassian filter of width="+str(z_sigma)
            data1 = scipy.ndimage.filters.gaussian_filter1d(data1, sigma=z_sigma, axis=2)

        
    fig = plt.figure(dpi=72)
    sub_fig = fig.add_subplot(111)

    # extract a slice from the 3D cube
    x_index = DIM/2
    if args.zindex < 0:
        print "Taking a yz slice at x index="+str(x_index)
        slice = data1[x_index,:,:]
        endstr = '_xindex'+str(x_index)
    else:
        print "Taking an xy slice at z index="+str(z_index)
        slice = data1[:,:,z_index]
        endstr = '_zindex'+str(z_index)
        if args.delzindex >= 0: #difference image is wanted
            del_z_index = args.delzindex
            other_z_index = int(z_index+del_z_index)
            print "Subtracting the slice at index="+str(other_z_index)
            slice = slice - data1[:,:,other_z_index]
            endstr = '_zindex'+str(z_index)+'-'+str(z_index+del_z_index)



    
    # check box type to determine default plotting options
    # check if it is a 21cm brightness temperature box
    if basename(filename)[0:3]=='del':
        if args.min > 1e4:
            minrange = -120
        if args.max < -1e4:
            maxrange = 20
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
        if args.min > 1e4:
            minrange = 0
        if args.max < -1e4:
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
        if args.min > 1e4:
            minrange = -0.5
        if args.max < -1e4:
            maxrange = 0.5
        slice = np.log10(1+slice)
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

        # check it is a recombination box
    elif basename(filename)[0:3]=='Nre':
        if args.min > 1e4:
            minrange = 0.0
        if args.max < -1e4:
            maxrange = 3
        cmap = LinearSegmentedColormap.from_list('mycmap', ['black', 'darkblue', 'yellow', 'red'])
        norm = MidpointNormalize(midpoint=1.5)
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
        c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm N_{rec}}$', fontsize=24, rotation=-90, labelpad=32)
        tick_array = np.linspace(minrange, maxrange, 4)

       # check it is a Gamma12 box
    elif basename(filename)[0:3]=='Gam':
        if args.min > 1e4:
            minrange = -2
        if args.max < -1e4:
            maxrange = 0
        slice = np.log10(slice)
        cmap = LinearSegmentedColormap.from_list('mycmap', ['black', 'darkblue', 'yellow', 'red'])
        norm = MidpointNormalize(midpoint=-1)
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
        c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
        c_dens.set_clim(vmin=minrange,vmax=maxrange)
        c_bar = fig.colorbar(c_dens, orientation='vertical')
        c_bar.set_label(r'${\rm log_{10}(\Gamma_{12}) }$', fontsize=24, rotation=-90, labelpad=32)
        tick_array = np.linspace(minrange, maxrange, 5)

    c_bar.set_ticks(tick_array)
	
    for t in c_bar.ax.get_yticklabels():
        t.set_fontsize(14)

    plt.savefig(filename+endstr+'.pdf', bbox_inches='tight')

    # do we want to save the array file?
    if args.savearray:
        np.save(filename+endstr, slice)

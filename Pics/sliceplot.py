#!/usr/bin/env python
#Written by- Arpan Das
#Usage> python sliceplot.py 'filename'
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
#import glob
import os
import sys

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

for i in range(1,len(sys.argv)):
	path = sys.argv[i]
	filename="" + path.split("/")[-1]

	
	if basename(filename)[0:3]=='del':

		DIM = int("" + path.split("_")[-2])
		label=str("" + path.split("_")[-1])
		tick_array = np.linspace(-3, 0, 4)

		# read in the data cube located in 21cmFast/Boxes/delta_T*
		data1 = load_binary_data(path)

		data1.shape = (DIM, DIM, DIM)
		data1 = data1.reshape((DIM, DIM, DIM), order='F')


		fig = plt.figure(dpi=72)

		#plot a slice
		slice = data1[:,:,DIM/2]
		scipy.ndimage.filters.gaussian_filter(slice, sigma=1.0, output=slice)
		#slice = log10(slice)

		sub_fig = fig.add_subplot(111)
		#Change the color here if you want
		cmap = LinearSegmentedColormap.from_list('mycmap', ['yellow','red','black','green','blue'])
		norm = MidpointNormalize(midpoint=0)
		frame1 = plt.gca()
		frame1.axes.get_xaxis().set_ticks([])
		frame1.axes.get_yaxis().set_ticks([])
		frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
		c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
		c_dens.set_clim(vmin=-200,vmax=50)
		c_bar = fig.colorbar(c_dens, orientation='vertical')
		c_bar.set_label(r'${\rm \delta T_b [\mathrm{mK}]}$', fontsize=24, rotation=-90, labelpad=32)
		tick_array = np.linspace(-200,50, 6)
		c_bar.set_ticks(tick_array)
		#plt.suptitle('Delta T slice for $\\alpha = 0.8\, - \,T_{vir} = 10^4 \,K$', fontsize=15)
	
		for t in c_bar.ax.get_yticklabels():
     			t.set_fontsize(14)
		#plt.savefig('%s.png'%(filename), bbox_inches='tight')
		plt.savefig('%s.pdf'%(filename), bbox_inches='tight')

	elif basename(filename)[0:3]=='xH_':
		DIM = int("" + path.split("_")[-2])
		label=str("" + path.split("_")[-1])
		tick_array = np.linspace(-3, 0, 4)

		# read in the data cube located in 21cmFast/Boxes/delta_T*
		data1 = load_binary_data(path)

		data1.shape = (DIM, DIM, DIM)
		data1 = data1.reshape((DIM, DIM, DIM), order='F')


		fig = plt.figure(dpi=72)

		#plot a slice
		slice = data1[:,:,DIM/2]
		scipy.ndimage.filters.gaussian_filter(slice, sigma=1.0, output=slice)
		#slice = log10(slice)

		sub_fig = fig.add_subplot(111)
		#Change the color here if you want
		cmap = LinearSegmentedColormap.from_list('mycmap', ['white','black'])
		norm = MidpointNormalize(midpoint=0.5)
		frame1 = plt.gca()
		frame1.axes.get_xaxis().set_ticks([])
		frame1.axes.get_yaxis().set_ticks([])
		frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
		c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
		c_dens.set_clim(vmin=0,vmax=1.0)
		c_bar = fig.colorbar(c_dens, orientation='vertical')
		c_bar.set_label(r'${\rm x_{HI}}$', fontsize=24, rotation=-90, labelpad=32)
		tick_array = np.linspace(0,1.0, 6)
		c_bar.set_ticks(tick_array)
	
		for t in c_bar.ax.get_yticklabels():
     			t.set_fontsize(14)

		#plt.savefig('%s.png'%(filename), bbox_inches='tight')
		plt.savefig('%s.pdf'%(filename), bbox_inches='tight')

	elif basename(filename)[0:3]=='upd':
		
		DIM = int("" + path.split("_")[-2])
		label=str("" + path.split("_")[-1])
		tick_array = np.linspace(-3, 0, 4)
		# read in the data cube located in 21cmFast/Boxes/delta_T*
		data1 = load_binary_data(path)

		data1.shape = (DIM, DIM, DIM)
		data1 = data1.reshape((DIM, DIM, DIM), order='F')

		fig = plt.figure(dpi=72)

		#plot a slice
		slice = np.log10(1+data1[:,:,DIM/2])
		
		scipy.ndimage.filters.gaussian_filter(slice, sigma=1.0, output=slice)

		sub_fig = fig.add_subplot(111)
		#Change the color here if you want
		cmap = LinearSegmentedColormap.from_list('mycmap', ['darkblue', 'green', 'yellow', 'red'])
		norm = MidpointNormalize(midpoint=0)
		frame1 = plt.gca()
		frame1.axes.get_xaxis().set_ticks([])
		frame1.axes.get_yaxis().set_ticks([])
		frame1.set_xlabel(r'${\rm\longleftarrow %s \longrightarrow}$'%(label), fontsize=20)
		c_dens = sub_fig.imshow(slice,cmap=cmap,norm=norm)
		c_dens.set_clim(vmin=-0.5,vmax=0.5)
		c_bar = fig.colorbar(c_dens, orientation='vertical')
		c_bar.set_label(r'${\rm log(\Delta)}$', fontsize=24, rotation=-90, labelpad=32)
		tick_array = np.linspace(-0.5, 0.5, 5)
		c_bar.set_ticks(tick_array)
	
		for t in c_bar.ax.get_yticklabels():
     			t.set_fontsize(14)

		#plt.savefig('%s.png'%(filename), bbox_inches='tight')
		plt.savefig('%s.pdf'%(filename), bbox_inches='tight')



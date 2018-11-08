# Script to generate a summary plot which are lightcone slice, the averaged birightness temperature and the power spectrum at k = 0.1 Mpc^-1.
# USAGE: python LightConeSnapshot.py
# Before use this script, must generates the co-eval redshifts output by 21cmFAST/21CMMC.
# Simply run 'drive_logZscroll_Ts.c' first, then generate all the co-eval boxes and power spectra that are needed to run.
# NOTE: 
#     (1) This Script is optimized for the simulation box of 300 Mpc on a side with a 200^3 grid.
#     (2) If folder/file names are changed, specify the folder/file names, i.e. 'DataLocation_Fiducial' and 'DataLocation' etc.
import numpy
import math
import os
from scipy import interpolate
from decimal import *
import string
import pickle
import pylab as P
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage as ndimage
import glob
import re
from matplotlib.widgets import Slider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

C = 29979245800.0 # speed of light  (cm/s)
OMm = 0.31
OMl = 1.0 - OMm
hlittle = 0.68 # little hubble h
Ho = hlittle*3.2407e-18 # s^-1 at z=0

# function DTDZ returns the value of dt/dz at the redshift parameter z.
def dtdz(z):
    x = numpy.sqrt( OMl/OMm ) * pow(1+z, -3.0/2.0)
    dxdz = numpy.sqrt( OMl/OMm ) * pow(1+z, -5.0/2.0) * (-3.0/2.0)
    const1 = 2 * numpy.sqrt( 1 + OMm/OMl ) / (3.0 * Ho)
    
    numer = dxdz * (1 + x*pow( pow(x,2) + 1, -0.5))
    denom = x + numpy.sqrt(pow(x,2) + 1)
    
    return (const1 * numer / denom)
# comoving distance (in cm) per unit redshift
def drdz(z):
    return (1.0+z)*C*dtdz(z)

if __name__ == '__main__':

    # Generates the co-eval redshifts output by 21cmFAST/21CMMC.
	# Set z_min, z_max and z_step_factor. These must be the same as in 21cmFAST/21CMMC.
    ZLOW = 6. # ZLOW in Programs/drive_logZscroll_Ts.c
    ZHIGH = 30. # Z_HEAT_MAX in Parameter_files/HEAR_PARAMS.H
    z_step_factor = 1.02 # ZPRIME_STEP_FACTOR in Parameter_files/HEAT_PARAMS.H
    z = ZLOW*1.0001

    # Box size and resolution of the 21cm Boxes
    BoxRes = '200' # Sould be STRING
    BoxSize_Float = 300. # Sould be INTEGER
    BoxSize = str(int(BoxSize_Float)) + 'Mpc'

    Redshifts = []
    Redshifts_float = []
    count = 0
    while z < ZHIGH:
        count = count + 1
        zp = str(round(z,6))
        Redshifts.append(zp)
        Redshifts_float.append(round(z,6))
        z = ((1.+z)*z_step_factor - 1.)
    N_REDSHIFT = count
    Redshifts_float = Redshifts_float[::-1]
	
    # Generate redshifts corresponding to each pixel.
    dR = (BoxSize_Float / float(BoxRes)) * 3.086e24 # size of cell (in comoving cm), 3.086e24 -> cm/Mpc.
    z1_LC = Redshifts_float[N_REDSHIFT-1]
    z_LC = start_z = z1_LC
    i = 0
    LightConeRedshifts = [] 
    while z1_LC < Redshifts_float[0]:
        z2_LC = float(Redshifts_float[N_REDSHIFT-2-i])
        while z_LC < z2_LC: # until we move to the next set of boxes
            LightConeRedshifts.append(z_LC)
            z_LC -= dR / drdz(z_LC)
        z1_LC = z2_LC
        i += 1

    # Find the redshifts for the light-cone boxes. Specifying the start and end redshift of each component of the light-cone
    FolderName = "../Boxes"
    filename_full = glob.glob('%s/delta_T_v3__zstart*_zend*_FLIPBOXES0_%s_%s_lighttravel'%(FolderName,BoxRes,BoxSize))
    Redshifts_LightCone_Begin =[]
    Redshifts_LightCone_End = []
    for i in range(len(filename_full)):
        filename = os.path.basename(filename_full[i])
        Redshifts_LightCone_Begin.append(str(filename[18:-48]))
        Redshifts_LightCone_End.append(str(filename[32:-34]))
    Redshifts_LightCone_Begin.sort()
    Redshifts_LightCone_End.sort()
		
    # The total number of light-cone boxes
    NumBoxes = len(Redshifts_LightCone_Begin)

    # Beginning and end redshifts for the figures. Nothing interesting happens beyond z = 25, so I cut it there.
    zmin = ZLOW
    zmax = 25. # Check the z range of the hightest lightcone box is enough to capture this range. 
               # if the box size is quite big (e.g. 1.5Gpc) and lightcone boxes generated with 'ZHIGH = 30', this quantity should be reduced.

    # Determine the corresponding z pixel index with the redshift maximum set above
    for jj in range(len(LightConeRedshifts)):
        if numpy.fabs(LightConeRedshifts[jj] - zmax) <= 0.01:
            Nmax = jj + 1

    # A randomly chosen number (between 0 and the pixel size of the 21cm boxes). Since we can only plot a slice, must pick a slice number
    FixedYPos_coord = int(BoxSize_Float/2) # I just set an half of box size. One can set any pixel number.

    NUM_Z_PS = len(Redshifts)

    # Arbitrarily chosen min and max values for the 21cm PS and the global signal
    PSmin = 0.01
    PSmax = 1000.

    min_val = -150.
    max_val = 30.

    # Numpy arrays to hold the light-cone slice data
    LightConeBox = numpy.zeros((int(BoxRes),Nmax))
    LightConeBox_interp = numpy.zeros((int(BoxRes),Nmax))

    ###########################################################################################
    # Fiducial model
    ###########################################################################################
    """
    If want to compare the variation of each parameter relative to a fiducial model,
    provide the directory location of the global averaged data.
    All word ending with "_Fiducial" represent data for the fiducial model, 
	e.g. 'DataLocation_Fiducial', 'NeutralFraction_Ave_Fiducial' and etc.
	This is appeared with grey lines.
    """
    DataLocation_Fiducial = "../Programs"
    Filename_Fiducial = "Power_k0.1"
    Fullname_Fiducial = os.path.join(DataLocation_Fiducial,Filename_Fiducial)

    # Read in the global averaged quantities (Neutral fraction and brightness temperature)
    Redshift_Fiducial = numpy.loadtxt(Fullname_Fiducial, usecols=(0,))
    NeutralFraction_Ave_Fiducial = numpy.loadtxt(Fullname_Fiducial, usecols=(1,))
    AveBrightness_Fiducial = numpy.loadtxt(Fullname_Fiducial, usecols=(2,))

    # Provide the directory location of the 21cm PS data
    DataLocation_Fiducial = '../Output_files/Deldel_T_power_spec'

    # Read in the k-values
    k_values_Fiducial = numpy.loadtxt('%s/ps_z%06.2f_nf%.6f_useTs1_aveTb%06.2f_%s_%s_v3'%(DataLocation_Fiducial,float(Redshifts[0]),NeutralFraction_Ave_Fiducial[0],AveBrightness_Fiducial[0],BoxRes,BoxSize), usecols=(0,))

    # Define a numpy array to hold all the 21cm PS data
    PS_values_Fiducial = numpy.zeros((NUM_Z_PS,len(k_values_Fiducial))) 

    # Read in the 21cm PS data
    for ii in range(NUM_Z_PS):
        PS_values_Fiducial[ii] = numpy.loadtxt('%s/ps_z%06.2f_nf%.6f_useTs1_aveTb%06.2f_%s_%s_v3'%(DataLocation_Fiducial,float(Redshifts[ii]),NeutralFraction_Ave_Fiducial[ii],AveBrightness_Fiducial[ii],BoxRes,BoxSize), usecols=(1,))


    # Define another array which will be used for defining what function to interpolate
    PS_Spline_Fiducial = numpy.zeros(len(k_values_Fiducial))

    # Define some matplotlib global variables
    matplotlib.rcParams.update({'font.size': 12})
    matplotlib.rc('text',usetex=True)
    matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

    # Define the colour map to be used to present the 21cm light-cone data
    EoR_colour = matplotlib.colors.LinearSegmentedColormap.from_list('mycmap', [(0, 'white'),(0.33, 'yellow'),(0.5, 'orange'),(0.68, 'red'),(0.833, 'black'),(0.87, 'blue'),(1, 'cyan')])
    plt.register_cmap(cmap=EoR_colour)


    sub_adj_l = 0.06
    sub_adj_b = 0.07
    sub_adj_r = 0.95
    sub_adj_t = 0.97
    sub_adj_w = 0.05
    sub_adj_h = 0.05


    # Numpy arrays for splined quantities
    SplinedFunction = numpy.zeros(Nmax)
    SplinedFunction_Tb_Fiducial = numpy.zeros(Nmax)
    SplinedFunction_PS_Fiducial = numpy.zeros(Nmax)

    # Numpy arrays for the redshift/frequency axes of the figure
    LightConeRedshifts_ForFigure = numpy.zeros(Nmax)
    Frequencies_ForFigure = numpy.zeros(Nmax)

    # Numpy array for the PS vs z plot (evaluated at a single k)
    PS_Evolution = numpy.zeros(NUM_Z_PS)

    # Numpy arrays for the x and y axis
    Y_Vals = numpy.zeros(int(BoxRes))
    X_Vals = numpy.zeros(Nmax)

    # Box width for the y-axis
    for ii in range(int(BoxRes)):
        Y_Vals[ii] = (0. + BoxSize_Float*(float(ii) + 0.)/((int(BoxRes)-1)))

    # Calculating the x-axis quantities. Either redshift, frequency or pixel number
    for ii in range(Nmax):
        LightConeRedshifts_ForFigure[ii] = LightConeRedshifts[ii]
        Frequencies_ForFigure[ii] = ((2.99792*10**8)/(0.2112*(1. + LightConeRedshifts_ForFigure[ii])))/(10**6)
        X_Vals[ii] = ii

    # Set up the spline for the global signal
    splined_Fiducial = interpolate.splrep(Redshifts,AveBrightness_Fiducial,s=0)
    
    # Evaluate the spline for the global signal
    for jj in range(Nmax):
        SplinedFunction_Tb_Fiducial[jj] = interpolate.splev(LightConeRedshifts_ForFigure[jj],splined_Fiducial,der=0)

    # k-value for which we present the 21cm PS as a function of redshift
    k_val = 0.1

    # Loop through all co-eval boxes, and spline interpolate the corresponding 21cm PS at the chosen k-value
    # Prior to performing the spline, check to see if there are any 'nan' values in the 21cm PS. These can arise when the neutral fraction is zero
    # If any 'nan' values are encountered, just set the value to be very small
    for ii in range(NUM_Z_PS):

        for jj in range(len(k_values_Fiducial)):
            if numpy.isnan(PS_values_Fiducial[ii][jj]) == True:
                PS_Spline_Fiducial[jj] = 0.000001
            else:
                PS_Spline_Fiducial[jj] = PS_values_Fiducial[ii][jj]

        splined_Fiducial = interpolate.splrep(k_values_Fiducial,numpy.log10(PS_Spline_Fiducial),s=0)

        PS_Evolution[ii] = 10**(interpolate.splev(k_val,splined_Fiducial,der=0))    

    # Spline the corresponding PS vs z data (at fixed k)
    splined_Fiducial = interpolate.splrep(Redshifts,numpy.log10(PS_Evolution),s=0)
    
    for jj in range(Nmax):
        SplinedFunction_PS_Fiducial[jj] = 10**(interpolate.splev(LightConeRedshifts_ForFigure[jj],splined_Fiducial,der=0))


    ###########################################################################################
    # To compare with default model
    ###########################################################################################
    """
    If want to compare the variation of each parameter relative to a fiducial model,
    provide the directory location of the global averaged data.
	In default setting, this part is the same with fiducial model.
    """
    DataLocation = "../Boxes"
    for k in range(NumBoxes):

        # Create the file string
        BoxNames = "%s/delta_T_v3__zstart%s_zend%s_FLIPBOXES0_%s_%s_lighttravel"%(DataLocation,Redshifts_LightCone_Begin[k],Redshifts_LightCone_End[k],BoxRes,BoxSize)

        # Read in the 21cm brightness temperature data from file
        f = open("%s"%(BoxNames),'rb')
        IndividualLightConeBox = numpy.fromfile(f, dtype = numpy.dtype('float32'), count = int(BoxRes)*int(BoxRes)*int(BoxRes))    
        f.close()

        # Take out the randomly determined slice and store the full light-cone slice
        for ii in range(int(BoxRes)):
            for kk in range(int(BoxRes)):
                if ii+int(BoxRes)*k < Nmax:
                    LightConeBox[kk][ii+int(BoxRes)*k] = IndividualLightConeBox[ii + int(BoxRes)*( FixedYPos_coord + int(BoxRes)*kk )]


    # Directory of the global averaged data
    DataLocation = "../Programs"

    # Read in the global averaged data
    NeutralFraction_Ave = numpy.loadtxt('%s/Power_k0.1'%(DataLocation), usecols=(1,))
    AveBrightness = numpy.loadtxt('%s/Power_k0.1'%(DataLocation),usecols=(2,))

    # Directory of the 21cm PS data
    DataLocation = "../Output_files/Deldel_T_power_spec"

    # Read in the k-values of the 21cm PS
    k_values = numpy.loadtxt('%s/ps_z%06.2f_nf%.6f_useTs1_aveTb%06.2f_%s_%s_v3'%(DataLocation_Fiducial,float(Redshifts[0]),NeutralFraction_Ave_Fiducial[0],AveBrightness_Fiducial[0],BoxRes,BoxSize), usecols=(0,))

    # Create the numpy array to contain all the 21cm PS data
    PS_values = numpy.zeros((NUM_Z_PS,len(k_values))) 

    # Read in the 21cm PS data
    for ii in range(NUM_Z_PS): # This is correct one
        #PS_values[ii] = numpy.loadtxt('%s/ps_z%06.2f_nf%.6f_useTs1_aveTb%06.2f_%s_%s_v3'%(DataLocation,float(Redshifts[ii]),NeutralFraction_Ave[ii],AveBrightness[ii],BoxRes,BoxSize), usecols=(1,))
        PS_values[ii] = numpy.loadtxt('%s/ps_z%06.2f_nf%.6f_useTs1_aveTb%06.2f_%s_%s_v3'%(DataLocation,float(Redshifts[ii]),NeutralFraction_Ave[ii],AveBrightness[ii],BoxRes,BoxSize), usecols=(1,))


    # Create a numpy array for spline interpolation of the 21cm PS
    PS_Spline = numpy.zeros(len(k_values)) 

    ###################################
    # Create the figure environment
    ###################################

    # Setup quantities for the plotting of the information.
    # Size of the y-panel
    ypanel_min = 0.0
    ypanel_max = 1.0

    # Background colour for the new panels to be created
    axcolor = 'white'

    # ********** set the panel position *************************************
    # - First, set margins    
    y_margin_top = 0.12 # margin at the top of the panel
    y_margin_btm = 0.12 # margin at the bottom of the panel
    x_margin_lft = 0.054 # margin at the left side of the panel
    x_margin_rgt = 0.05 # margin at the right side of the panel

    N_main_panels = 3 # the number of main panels which are lightcone, the brightness temperature and PS at k=0.1 Mpc^-1.

    # - Then, set the width and height of panels.
    ratio_main_panel = 1.0 # want use 100% of overall width for main panels
    panel_width = (1.0 - x_margin_lft - x_margin_rgt)*ratio_main_panel
    panel_height = (1.0 - y_margin_top - y_margin_btm)/float(N_main_panels)
    panel_xpos = x_margin_lft

    panel_ypos =  1.0 - y_margin_top - panel_height 

    panel2_xpos = panel_xpos    
    panel2_ypos = panel_ypos - panel_height 
    panel2_width = panel_width
    panel2_height = panel_height

    panel3_xpos = panel_xpos
    panel3_ypos = panel2_ypos - panel_height 
    panel3_width = panel_width
    panel3_height = panel_height

    # To avoid distortion of light-cone panel, 
    # the width of the panel has to be the same with the height of the panel times the number of Light-cone boxes (NumBoxes), 
    # i.e. NumBoxes * panel_height * overall y length = panel_width * overall x length
    x_total = 18. # set an arbitrary width of the figure
    y_total = panel_width*x_total/float(NumBoxes)/panel_height
    print 'figure size is (',x_total,',',y_total,')'
    fig = plt.figure(figsize=(x_total,y_total))


###############################################################################################################################
#       For creating the light-cone panel
###############################################################################################################################
    ax = plt.axes([panel_xpos, panel_ypos, panel_width, panel_height])

    plt.subplots_adjust(sub_adj_l,sub_adj_b,sub_adj_r,sub_adj_t,sub_adj_w,sub_adj_h)
    # Add the x and y labels
    ax.set_ylabel(r"$\boldsymbol{L\,{\rm [Mpc]}}$",fontsize='large',labelpad=15)
    ax.set_xlabel(r"$\boldsymbol{{\rm Frequency\,[MHz]}}$",fontsize='large',labelpad=10)
    ax.xaxis.set_label_position('top') 

    ax.axis([X_Vals[0],X_Vals[-1],Y_Vals[0],Y_Vals[-1]])
    
    # Plot the light-cone data
    CS = ax.pcolormesh(X_Vals,Y_Vals,LightConeBox, vmin=min_val, vmax=max_val,cmap=EoR_colour,shading='gouraud')

    # Create the colour bar, and normalise the colour bar to be between min_val and min_val
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="1%", pad=0.05, pack_start=False)
    fig.add_axes(cax)

    m = plt.cm.ScalarMappable(cmap=EoR_colour)
    m.set_array(LightConeBox)
    m.set_clim(vmin=min_val, vmax=max_val)

    # Add the colour bar, and the label
    cbar = plt.colorbar(m,cax=cax)
    cbar.set_label(r"$\boldsymbol{\delta T_{\rm b}\,[{\rm mK}]}$",rotation=270,labelpad=8)

    # Place the labelling and ticks for the colour bar
    #cbar.set_ticks(numpy.arange(-150., 31., 30.))
    cbar.set_ticks(numpy.arange(-150., 31., 30.))
    cbar.set_ticklabels([r"$-150$",r"$-120$",r"$-90$",r"$-60$",r"$-30$",r"$0$",r"$30$"])
    cbar_yticks = cbar.ax.yaxis.get_ticklabels()
    cbar_yticks[1].set_visible(False)
    cbar_yticks[3].set_visible(False)

    cbar.ax.tick_params(axis='y', direction='out',pad=2)


    # Manually add the tick locations for the frequency axis
    frequency_locations_z = numpy.arange(60.,201.,10.)
    frequency_locations_minor_z = numpy.arange(55.,205.,10.)
    frequency_locations_z = frequency_locations_z[::-1]
    frequency_locations_minor_z = frequency_locations_minor_z[::-1]

    # Add the labels for the frequency axis
    frequency_tick_labels = [r"$200$",r"$190$",r"$180$",r"$170$",r"$160$",r"$150$",r"$140$",r"$130$",r"$120$",r"$110$",r"$100$",r"$90$",r"$80$",r"$70$",r"$60$"]

    frequency_tick_positions_z = []

    # Find the corresponding pixel index to the frequency tick locations (major) to be added.
    for ii in range(len(frequency_locations_z)):

        for jj in range(len(Frequencies_ForFigure)):
            if Frequencies_ForFigure[jj] <= frequency_locations_z[ii]:
                frequency_tick_positions_z.append(jj)
                break

    # Find the corresponding pixel index to the frequency tick locations (major) to be added.
    frequency_tick_positions_minor_z = []

    for ii in range(len(frequency_locations_minor_z)):

        for jj in range(len(Frequencies_ForFigure)):
            if Frequencies_ForFigure[jj] <= frequency_locations_minor_z[ii]:
                frequency_tick_positions_minor_z.append(jj)
                break

    # Setup the x-axis ticks, corresponding to the values and labels above
    ax.set_xlim(0, Nmax)
    ax.set_xticks(frequency_tick_positions_z)
    ax.set_xticklabels(frequency_tick_labels)
    ax.set_xticks(frequency_tick_positions_minor_z,minor=True)


    # Set up the y-axis tick locations and labels
    ytick_locations = numpy.arange(0., 301., 50.)
    ytick_labels = [r"$0$",r"$50$",r"$100$",r"$150$",r"$200$",r"$250$",r"$300$"]

    ax.set_yticks(ytick_locations)
    ax.set_yticklabels(ytick_labels)
    yticks = ax.set_yticklabels(ytick_labels)
    yticks[0].set_visible(False)
    yticks[1].set_visible(False)
    yticks[3].set_visible(False)
    yticks[5].set_visible(False)

    # Set the y-axis minor tick labels
    minorLocator = matplotlib.ticker.MultipleLocator(10.)
    ax.yaxis.set_minor_locator(minorLocator)

    # Turn of the x-axis labels at the bottom and turn on the x-axis labels at the top
    ax.tick_params(axis='x',direction='out', which='both', bottom='off', top='on',  labelbottom='off',  labeltop='on', pad=2)
    ax.tick_params(axis='y',direction='out', which='both', left='on', right='off',  labelleft='on',pad=2)

    # Increase the size, change the text stype and bold the text of the tick marks (for the colour bar)
    for tick in cbar.ax.yaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

    # Increase the size, change the text stype and bold the text of the tick marks (for the x-axis)
    for tick in ax.xaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

    # Increase the size, change the text stype and bold the text of the tick marks (for the y-axis)
    for tick in ax.yaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

###############################################################################################################################
#       For creating the global signal panel
###############################################################################################################################
    # Panel location
    ax = plt.axes([panel2_xpos, panel2_ypos, panel2_width, panel2_height])
    ax.set_ylabel(r"$\boldsymbol{\overline{\delta T}_{\rm b}\,\,[\rm mK]}$",fontsize='large',labelpad=0)

    plt.subplots_adjust(sub_adj_l,sub_adj_b,sub_adj_r,sub_adj_t,sub_adj_w,sub_adj_h)

    # Range for panel
    #ax.axis([0,Nmax,min_val,max_val])
    ax.axis([0,Nmax,-120,20])


    # Prepere the spline for the global signal
    splined = interpolate.splrep(Redshifts,AveBrightness,s=0)

    for jj in range(Nmax):
        SplinedFunction[jj] = interpolate.splev(LightConeRedshifts_ForFigure[jj],splined,der=0)

    # Plot the global signal value of the fiducial setup (i.e. this is always fixed)
    lines = plt.plot(X_Vals,SplinedFunction_Tb_Fiducial,'-')
    plt.setp(lines, color='grey', linewidth=2.0)
 
    # Plot the global signal for this specific varied parameter
    lines = plt.plot(X_Vals,SplinedFunction,'-')
    plt.setp(lines, color='black', linewidth=3.0)
 
    # Horizontal line highlighting the zero line
    ax.plot((X_Vals[0],X_Vals[-1]),(0.,0.),color='black',ls='dashed',linewidth=1.)

    # Add a colour bar, but only to ensure that the width of the figure panel is consistent with the light-cone strip above it
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="1%", pad=0.05, pack_start=False)
    fig.add_axes(cax,frameon=False)
    cax.axis('off')

    cax.tick_params(axis='y',direction='out', which='both', left='off', right='off',  labelleft='off')
    cax.tick_params(axis='x',direction='out', which='both', top='off', bottom='off',  labelbottom='off')

    # Manually add the major tick locations for the redshift axis.
    xtick_locations_z = numpy.arange(6.,25.1,1)
    xtick_labels_z = [r"$6.0$",r"$7.0$",r"$8.0$",r"$9.0$",r"$10.0$",r"$\,$",r"$12.0$",r"$\,$",r"$14.0$",r"$\,$",r"$16.0$",r"$\,$",r"$18.0$",r"$\,$",r"$20.0$",r"$\,$",r"$\,$",r"$\,$",r"$\,$",r"$25.0$"]
    #xtick_labels_z = [r"$6.0$",r"$7.0$",r"$8.0$",r"$9.0$",r"$10.0$",r"$11.0$",r"$12.0$",r"$13.0$",r"$14.0$",r"$15.0$",r"$16.0$",r"$17.0$",r"$18.0$",r"$19.0$",r"$20.0$",r"$21.0$",r"$22.0$",r"$23.0$",r"$24.0$",r"$25.0$"]

    xtick_positions_z = []

    for ii in range(len(xtick_locations_z)):

        for jj in range(len(LightConeRedshifts_ForFigure)):
            if LightConeRedshifts_ForFigure[jj] >= xtick_locations_z[ii]:
                xtick_positions_z.append(jj)
                break

    ax.set_xlim(0, Nmax)
    ax.set_xticks(xtick_positions_z)
    ax.set_xticklabels(xtick_labels_z)
    xticks = ax.set_xticklabels(xtick_labels_z)

    # Add the minor tick locations
    xtick_minor_locations_z = numpy.arange(6.5,24.6,1.)

    xtick_minor_positions_z = []

    for ii in range(len(xtick_minor_locations_z)):

        for jj in range(len(LightConeRedshifts_ForFigure)):
            if LightConeRedshifts_ForFigure[jj] >= xtick_minor_locations_z[ii]:
                xtick_minor_positions_z.append(jj)
                break

    ax.set_xticks(xtick_minor_positions_z,minor=True)

    # Add the minor tick locations for the y-axis
    minorLocator = matplotlib.ticker.MultipleLocator(10.)
    ax.yaxis.set_minor_locator(minorLocator)

    ytick_locations = numpy.arange(-120.,21.,20.)
    ytick_labels = [r"$-120$",r"$-100$",r"$-80$",r"$-60$",r"$-40$",r"$-20$",r"$0$",r"$20$"]

    ax.set_yticks(ytick_locations)
    yticks = ax.set_yticklabels(ytick_labels)
    yticks[1].set_visible(False)
    yticks[3].set_visible(False)
    yticks[5].set_visible(False)
    yticks[7].set_visible(False)

    ax.tick_params(axis='y',direction='out', which='both', left='on', right='on',  labelleft='on',pad=2)
    ax.tick_params(axis='x',direction='in', which='both', top='on', bottom='on',  labelbottom='off')

    for tick in ax.xaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

    for tick in ax.yaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

###############################################################################################################################
#       Now add the final panel, which is the 21cm PS as a function of redshift for a single k-value
###############################################################################################################################

    ax = plt.axes([panel3_xpos, panel3_ypos, panel3_width, panel3_height])
    # Add axis labels
    ax.set_xlabel(r"$\boldsymbol{{\rm Redshift},\,z}$",fontsize='large',labelpad=8)
    ax.set_ylabel(r"$\boldsymbol{\overline{\delta T}^{2}_{\rm b}\Delta^{2}_{21}\,\,[\rm mK^{2}]}$",fontsize='large',labelpad=1)

    ax.axis([zmin,zmax,PSmin,PSmax])

    ax.set_yscale('log')

    k_val = 0.1

    # Spline evaluate the 21cm PS at each redshift to determine the PS value at k_val
    # Prior to calculating the spline, check to see if any 'nan' values exists, and if so
    # replace them with a small value
    for ii in range(NUM_Z_PS):

        for jj in range(len(k_values)):
            if numpy.isnan(PS_values[ii][jj]) == True:
                PS_Spline[jj] = 0.000001
            elif (PS_values[ii][jj] < 0.000001):
                PS_Spline[jj] = 0.000001
            else:
                PS_Spline[jj] = PS_values[ii][jj]

        splined = interpolate.splrep(k_values,numpy.log10(PS_Spline),s=0)

        PS_Evolution[ii] = 10**(interpolate.splev(k_val,splined,der=0))

    # Spline the 21cm PS as a function of redshift
    splined = interpolate.splrep(Redshifts,numpy.log10(PS_Evolution),s=0) 

    for jj in range(Nmax): 
        SplinedFunction[jj] = 10**(interpolate.splev(LightConeRedshifts_ForFigure[jj],splined,der=0))


    # Plot the fiducial 21cm PS as a function of redshift
    lines = plt.plot(X_Vals,SplinedFunction_PS_Fiducial,'-')
    plt.setp(lines, color='grey', linewidth=2.0)

    # Plot the 21cm PS as a function of redshift for this parameter set
    plt.setp(lines, color='black', linewidth=3.0)

    # Manually place the major y-axis tick marks and corresponding labels
    ytick_locations = [0.01,0.1,1.0,10.,100.,1000.]
    ytick_labels = [r"$10^{-2}$",r"$10^{-1}$",r"$10^{0}$",r"$10^{1}$",r"$10^{2}$",r"$10^{3}$"]

    ax.set_yticks(ytick_locations)
    yticks = ax.set_yticklabels(ytick_labels)
    yticks[1].set_visible(False)
    yticks[3].set_visible(False)
    yticks[5].set_visible(False)

    # Place a text label corresponding to the k-value at which we have evaluated the 21cm PS
    plt.text(xtick_positions_z[15], 150.0, r"$\boldsymbol{k = 0.1\,{\rm Mpc}^{-1}}$",color='black', fontsize=12)


    # Add in the x-tick locations
    ax.set_xlim(0, Nmax)
    ax.set_xticks(xtick_positions_z)
    xticks = ax.set_xticklabels(xtick_labels_z)
    xticks[5].set_visible(False)
    xticks[7].set_visible(False)
    xticks[9].set_visible(False)
    xticks[11].set_visible(False)
    xticks[13].set_visible(False)
    xticks[15].set_visible(False)
    xticks[16].set_visible(False)
    xticks[17].set_visible(False)
    xticks[18].set_visible(False)

    ax.set_xticks(xtick_minor_positions_z,minor=True)

    ax.tick_params(axis='y',direction='out', which='both', left='on', right='on',  labelleft='on',pad=2)
    ax.tick_params(axis='x',direction='in', which='both', top='on', bottom='on',  labelbottom='on',pad=3)

    # Again, add a colour bar purely for ensuring the correct figure width
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size="1%", pad=0.05, pack_start=False)
    fig.add_axes(cax,frameon=False)
    cax.axis('off')

    cax.tick_params(axis='y',direction='out', which='both', left='off', right='off',  labelleft='off')
    cax.tick_params(axis='x',direction='out', which='both', top='off', bottom='off',  labelbottom='off')

    for tick in ax.xaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

    for tick in ax.yaxis.get_ticklabels():
        tick.set_fontsize('large')
        tick.set_fontname('Times New Roman')
        tick.set_weight('bold')

            

    # The file format is required to have the file name with an incrementally increasing index
    plt.savefig('LightConeSnapshot.png',dpi=1000)

    plt.close()

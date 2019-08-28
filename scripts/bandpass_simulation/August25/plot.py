# ----------------------------------------------------------------------------
#
# TITLE -
# AUTHOR - James Lane
# PROJECT -
# CONTENTS:
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Plot the results of the August 25 simulations
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb
# import copy
# import glob
# import subprocess

## Plotting
from matplotlib import pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import colors
# from matplotlib import cm

## scipy
import scipy.io as sio

# ----------------------------------------------------------------------------

plt.rc('text',usetex=True)

# ----------------------------------------------------------------------------

def load_matlab_array(filename,dict_key=None):
    '''load_matlab_array:

    Load data from a single .mat format file

    Args:

    Returns:
        data_arr
    '''
    data_in = sio.loadmat(filename)

    # If no key was supplied then search for a single non-trivial
    # dictionary key
    if dict_key == None:
        for key in data_in.keys():
            if dict_key != None: raise RuntimeError('Too many non-trivial keys!')
            if key[:2] != '__' and key[-2:] != '__':
                dict_key = key
            ##fi
        ###i
    ##fi

    data_out = data_in[dict_key]
    return data_out.reshape(data_out.shape[1])
#def

# ----------------------------------------------------------------------------

bands = ['zJ','YJ','JH','HKs']
mags = [6,8,10,12,14]
mods = [1,3,5]
band_colors = ['DodgerBlue','ForestGreen','DarkOrange','Red']
mod_linestyles = ['solid','dashed','dotted']

fig = plt.figure()
ax = fig.add_subplot(111)

# Loop over all bands and magnitudes
for i in range( len(bands) ):

    for j in range( len(mods) ):
        
        # Array to hold strehls for a given band and magnitude
        strehl_band = np.zeros( len(mags) )
        strehl_error = np.zeros( len(mags) )
        
        for k in range( len(mags) ):
            filename = './output/py_'+str(mags[k])+'mag_'+str(mods[j])+'mod_'+bands[i]+'band.mat'
            strehl = load_matlab_array(filename,dict_key='strehlVecMaster')
            strehl_band[k] = np.average(strehl)
            strehl_error[k] = np.std(strehl)
        ###k
        
        ax.errorbar( mags, 100*strehl_band, yerr=100*strehl_error, fmt='o',
            mec='Black', mfc=band_colors[i], ls=mod_linestyles[j], ms=1, 
            linewidth=1.0,
            c=band_colors[i], ecolor='Black', elinewidth=0.5, capsize=1, 
            capthick=0.5)
        
    ###j
###i

for i in range( len(bands) ):
    ax.plot([],[],color=band_colors[i],label=bands[i]+' Band')
###i
for i in range( len(mods) ):
    ax.plot([],[],color='Black',linestyle=mod_linestyles[i],
        label=str(mods[i])+r' $\lambda$/D')
###i

ax.set_ylim(50,95)
ax.set_xlim(9.5,14.5)
ax.set_xlabel('Magnitude')
ax.set_ylabel(r'$K-$band Strehl [\%]')
ax.legend()

plt.savefig('StrehlStats.pdf')
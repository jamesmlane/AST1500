# ----------------------------------------------------------------------------
#
# TITLE - plot.py
# AUTHOR - James Lane
# PROJECT - AST 1500
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Plot the results of the NGS bandpass simulations
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb, glob

## Plotting
from matplotlib import pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import colors
# from matplotlib import cm

## scipy
import scipy.io as sio

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

bands = ['J','H','K','JH','HK','JHK']
mags = [6,8,10,12,14]
band_colors = ['Orange','Coral','OrangeRed','LightBlue','DodgerBlue','Blue']

fig = plt.figure()
ax = fig.add_subplot(111)

# Loop over all bands and magnitudes
for i in range( len(bands) ):

    # Array to hold strehls for a given band
    strehl_band = np.zeros(len(mags))

    for j in range( len(mags) ):

        filename = './output/py_'+str(mags[j])+'mag_'+bands[i]+'band_scaledmodulation2.mat'
        strehl = load_matlab_array(filename,dict_key='strehlVecMaster')
        strehl_band[j] = np.average(strehl)

    ###j

    ax.plot( mags, 100*strehl_band, color=band_colors[i], label=bands[i]+'-band',
        marker='o', markeredgecolor='Black')
###i

ax.set_ylim(25,100)
ax.set_xlabel('Magnitude')
ax.set_ylabel('Strehl [%]')
ax.legend()

plt.savefig('StrehlStats.pdf')













# ----------------------------------------------------------------------------

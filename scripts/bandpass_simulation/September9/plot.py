# ----------------------------------------------------------------------------
#
# TITLE -
# AUTHOR - James Lane
# PROJECT -
# CONTENTS:
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Plot the results of the September 1 simulations
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb

## Plotting
from matplotlib import pyplot as plt

## scipy
sys.path.insert(0,'../../../src/')
import ast1500.io

# ----------------------------------------------------------------------------

plt.rc('text',usetex=True)

# ----------------------------------------------------------------------------

bands = ['z','Y','J','H','Ks']
mags = [10,12,14]
mods = [1,3,5]
band_colors = ['DodgerBlue','ForestGreen','DarkOrange','Red','Purple']
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
            filename = './output/py_'+str(mags[k])+'mag_'+str(mods[j])+'mod_03gain_'+bands[i]+'band.mat'
            strehl = ast1500.io.load_matlab_array(filename,
                dict_key='strehlVecMaster')
            strehl_band[k] = np.average(strehl)
            strehl_error[k] = np.std(strehl)
        ###k
        
        ax.errorbar( mags, 100*strehl_band, yerr=100*strehl_error, fmt='o',
            mec='Black', mfc=band_colors[i], ls=mod_linestyles[j], ms=4, 
            linewidth=1.0,
            c=band_colors[i], ecolor='Black', elinewidth=0.5, capsize=2, 
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

ax.set_ylim(40,95)
ax.set_xlim(9.5,14.5)
ax.set_xlabel('Guide Star Magnitude')
ax.set_ylabel(r'$K-$band Strehl [\%]')
ax.legend()

plt.savefig('StrehlStats_SingleBands.pdf')

# -------------

bands = ['zJ','YJ','JH','HKs']
mags = [10,12,14]
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
            filename = './output/py_'+str(mags[k])+'mag_'+str(mods[j])+'mod_03gain_'+bands[i]+'band.mat'
            strehl = ast1500.io.load_matlab_array(filename,
                dict_key='strehlVecMaster')
            strehl_band[k] = np.average(strehl)
            strehl_error[k] = np.std(strehl)
        ###k
        
        ax.errorbar( mags, 100*strehl_band, yerr=100*strehl_error, fmt='o',
            mec='Black', mfc=band_colors[i], ls=mod_linestyles[j], ms=4, 
            linewidth=1.0,
            c=band_colors[i], ecolor='Black', elinewidth=0.5, capsize=2, 
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

ax.set_ylim(40,95)
ax.set_xlim(9.5,14.5)
ax.set_xlabel('Guide Star Magnitude')
ax.set_ylabel(r'$K-$band Strehl [\%]')
ax.legend()

plt.savefig('StrehlStats_DoubleBands.pdf')
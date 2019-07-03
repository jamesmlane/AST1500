# ----------------------------------------------------------------------------
#
# TITLE - assess_ESO_spectra_normalization.py
# AUTHOR - James Lane
# PROJECT - AST1500
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Plot how well the normalization has done
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb

## Plotting
import matplotlib
from matplotlib import pyplot as plt

## Astropy
from astropy.io import fits

# Project specific
sys.path.append('../../src')
import ast1500.spectra

matplotlib.rc('text',usetex=True)

# ----------------------------------------------------------------------------

def calculate_magnitude(data,response,wavelength,dlambda):
    '''
    Assume dlambda is in angstrom and therefore c will be in angstroms
    '''
    
    c_angstrom = 3E18
    calc_num = np.sum( data*response*dlambda )
    calc_denom = np.sum( response*c_angstrom*dlambda/(wavelength**2) )    

    return -2.5*np.log10( calc_num / calc_denom )
#def

# ----------------------------------------------------------------------------

# First read the spectrum and the information file
spec_dir = '../../data/spectra/eso_spectral_lib/proper/'

spec_library = ast1500.spectra.ESOSpectralLibrary(path_to_spectra=spec_dir)
u_band_filter = ast1500.spectra.PhotometricFilter('U')
b_band_filter = ast1500.spectra.PhotometricFilter('B')
v_band_filter = ast1500.spectra.PhotometricFilter('V')
r_band_filter = ast1500.spectra.PhotometricFilter('R')
i_band_filter = ast1500.spectra.PhotometricFilter('I')
j_band_filter = ast1500.spectra.PhotometricFilter('J')
h_band_filter = ast1500.spectra.PhotometricFilter('H')
ks_band_filter = ast1500.spectra.PhotometricFilter('Ks')

spec_info_filename = '../../data/spectra/eso_spectral_lib/information.fit'
hdu = fits.open(spec_info_filename)
data = hdu[1].data

# Get information about the spectrum
Mv_sptype = data['SpType']
Mv = data['Mbol']-data['BCv']
Mi = data['Mbol']-data['BCi']
Mk = data['Mbol']-data['BCk']
vmj = data['V-J']
vmh = data['V-H']

fig = plt.figure(figsize=(5,15))
axs = fig.subplots(nrows=5,ncols=1)

# Loop over each spectrum
for i in range( spec_library.n_spectra ):
    
    # Read in the spectrum
    spec_name = spec_library.spectra_names[i][2:]
    spec_filename = spec_library.spectra_filenames[i]
    spec_wavelength, spec_data = np.genfromtxt(spec_filename,usecols=(0,1)).T
    
    # Figure out where this spectrum is in the information file
    where_spectra_in_info = np.where( Mv_sptype.upper() == spec_name.upper() )[0]
    this_Mv = Mv[where_spectra_in_info]
    this_Mi = Mi[where_spectra_in_info]
    this_Mk = Mk[where_spectra_in_info]
    this_vmj = vmj[where_spectra_in_info]
    this_vmh = vmh[where_spectra_in_info]
    
    # Calculate the response and difference
    v_band_response = v_band_filter.response(spec_wavelength/10.)
    i_band_response = i_band_filter.response(spec_wavelength/10.)
    j_band_response = j_band_filter.response(spec_wavelength/10.)
    h_band_response = h_band_filter.response(spec_wavelength/10.)
    ks_band_response = ks_band_filter.response(spec_wavelength/10.)
    delta_lambda = np.average( np.diff(spec_wavelength) )
    
    # Calculate the normalized V-band magnitude
    v_band_magnitude = calculate_magnitude(spec_data, v_band_response, 
        spec_wavelength, delta_lambda)
    i_band_magnitude = calculate_magnitude(spec_data, i_band_response, 
        spec_wavelength, delta_lambda)
    j_band_magnitude = calculate_magnitude(spec_data, j_band_response, 
        spec_wavelength, delta_lambda)
    h_band_magnitude = calculate_magnitude(spec_data, h_band_response, 
        spec_wavelength, delta_lambda)
    ks_band_magnitude = calculate_magnitude(spec_data, ks_band_response, 
        spec_wavelength, delta_lambda)
    
    axs[0].scatter( v_band_magnitude, this_Mv-v_band_magnitude, s=5, color='Black') 
    axs[1].scatter( i_band_magnitude, this_Mi-i_band_magnitude, s=5, color='Black') 
    axs[2].scatter( ks_band_magnitude, this_Mk-ks_band_magnitude, s=5, color='Black')
    axs[3].scatter( v_band_magnitude-j_band_magnitude,
        this_vmj - (v_band_magnitude-j_band_magnitude), s=5, color='Black')
    axs[4].scatter( v_band_magnitude-h_band_magnitude,
        this_vmh - (v_band_magnitude-h_band_magnitude), s=5, color='Black')
###i

# axs[0].plot([-10,10], [-10,10], linewidth=0.5, linestyle='dashed', color='Black')
# axs[1].plot([-10,10], [-10,10], linewidth=0.5, linestyle='dashed', color='Black')
# axs[2].plot([-10,10], [-10,10], linewidth=0.5, linestyle='dashed', color='Black')
# axs[0].axhline(0, linewidth=0.5, linestyle='dashed', color='Black')
axs[0].axhline(0.0, linewidth=0.5, linestyle='dashed', color='Red')
# axs[1].axhline(0, linewidth=0.5, linestyle='dashed', color='Black')
axs[1].axhline(-0.45, linewidth=0.5, linestyle='dashed', color='Red')
# axs[2].axhline(0, linewidth=0.5, linestyle='dashed', color='Black')
axs[2].axhline(-1.85, linewidth=0.5, linestyle='dashed', color='Red')
# axs[3].axhline(0, linewidth=0.5, linestyle='dashed', color='Black')
axs[3].axhline(0.91, linewidth=0.5, linestyle='dashed', color='Red')
# axs[4].axhline(0, linewidth=0.5, linestyle='dashed', color='Black')
axs[4].axhline(1.39, linewidth=0.5, linestyle='dashed', color='Red')

axs[0].set_xlim(-12,12)
# axs[0].set_ylim(-12,12)
axs[1].set_xlim(-12,12)
# axs[1].set_ylim(-12,12)
axs[2].set_xlim(-12,12)
# axs[2].set_ylim(-12,12)

axs[0].set_xlabel('Spectral $V$')
axs[0].set_ylabel('True $V$ $-$ Spectral $V$')
axs[1].set_xlabel('Spectral $I$')
axs[1].set_ylabel('True $I_{c}$ $-$ Spectral $I$')
axs[2].set_xlabel('Spectral $K_{s}$')
axs[2].set_ylabel('True $K$ $-$ Spectral $K_{s}$')
axs[3].set_xlabel('Spectral $(V-J)$')
axs[3].set_ylabel('True $(V-J)$ $-$ Spectral $(V-J)$')
axs[4].set_xlabel('Spectral $(V-H)$')
axs[4].set_ylabel('True $(V-H)$ $-$ Spectral $(V-H)$')

plt.tight_layout()
plt.savefig('magnitude_comparison.pdf')
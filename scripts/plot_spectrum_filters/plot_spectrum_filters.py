# ----------------------------------------------------------------------------
#
# TITLE - plot_spectrum_filters.py
# AUTHOR - James Lane
# PROJECT - AST1500
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Show some plots of spectra and filters
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

# First read the spectrum
spec_dir = '../../data/spectra/eso_spectral_lib/proper/'
spec_normalized = True
spec_library = ast1500.spectra.ESOSpectralLibrary(path_to_spectra=spec_dir)
spectral_names = ['ukg2v','ukk7v','ukm2v']
spectral_plot_names = [r'G2V',r'K7V',r'M2V']
n_spectra = len(spectral_names)
wavelength_range = [5000,25000] # angstroms
spec_color=['ForestGreen','DarkOrange','Red']

# ----------------------------------------------------------------------------

r_band_filter = ast1500.spectra.PhotometricFilter('R')
i_band_filter = ast1500.spectra.PhotometricFilter('I')
j_band_filter = ast1500.spectra.PhotometricFilter('J')
h_band_filter = ast1500.spectra.PhotometricFilter('H')
ks_band_filter = ast1500.spectra.PhotometricFilter('Ks')
bands = [r_band_filter,
         i_band_filter,
         j_band_filter,
         h_band_filter,
         ks_band_filter]
n_bands = len(bands)
band_plot_names=['R','I','J','H','Ks']
band_color=['Red','Orange','ForestGreen','DodgerBlue','Purple']

# ----------------------------------------------------------------------------

# Make the plot
fig = plt.figure(figsize=(8,8))
axs = fig.subplots(nrows=2,ncols=1)

for i in range(n_spectra):
    
    spec_wavelength, spec_data = spec_library.read_spectra(spectral_names[i],
        normalized=spec_normalized, return_data=True)
    where_in_wavelength_range = np.where( ( spec_wavelength > wavelength_range[0] ) &
                                          ( spec_wavelength < wavelength_range[1] ) )[0]
    spec_wavelength = spec_wavelength[where_in_wavelength_range]
    spec_data = spec_data[where_in_wavelength_range]
    spec_data /= np.max(spec_data)
    axs[0].plot(spec_wavelength, spec_data, color=spec_color[i], 
        label=spectral_plot_names[i])
###i
    
for j in range(n_bands):
    wavelength_arr = np.arange(wavelength_range[0],wavelength_range[1],1)
    band_response = bands[j].response(spec_wavelength)
    axs[1].plot(spec_wavelength, band_response, color=band_color[j],
        label=band_plot_names[j])
###j

axs[0].set_xlim(wavelength_range[0],wavelength_range[1])
axs[0].set_ylim(-0.1,1.1)
axs[1].set_xlim(wavelength_range[0],wavelength_range[1])
axs[1].set_ylim(0.01,1.1)

label_fs = 16
axs[0].set_xlabel(r'Wavelength [\AA]', fontsize=label_fs)
axs[0].set_ylabel(r'Flux [Normalized]', fontsize=label_fs)
axs[1].set_xlabel(r'Wavelength [\AA]', fontsize=label_fs)
axs[1].set_ylabel(r'Transmittance [Fraction]', fontsize=label_fs)

axs[0].legend()
axs[1].legend()

plt.savefig('spectrum_fiters_plot.pdf')
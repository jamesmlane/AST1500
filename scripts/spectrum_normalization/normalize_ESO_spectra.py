# ----------------------------------------------------------------------------
#
# TITLE - normalize_ESO_spectra.py
# AUTHOR - James Lane
# PROJECT - AST1500
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''Normalize ESO Pickles spectra so that they are expressed as Flambda and
return proper magnitudes. The default normalization is so that the spectra 
are equal to 1 at 5500 AA or something.
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb

## Plotting
from matplotlib import pyplot as plt

## Astropy
from astropy.io import fits

# Project specific
sys.path.append('../../src')
import ast1500.spectra

# ----------------------------------------------------------------------------

def calculate_magnitude(data,response,wavelength,dlambda):
    '''calculate_magnitude:
    
    Calculates the magnitude from flambda. Assume that flambda is in 
    erg / s / square cm / Angstrom and dlambda is in Angstroms. Therefore c 
    will be in Angstroms / s
    
    Only photometric system supported currently is AB magnitudes.
    
    If zero-point is None then a default will be used if possible.
    
    Args:
        response (float array)
        wavelength (float array)
        dlambda (float) - difference between wavelength bins
        photometric_system (str) - Name of photometric system ['AB']
        zero_point (float) - Zero point of the filter [None]
        
    Returns:
        magnitude (float)
    '''
    
    c_angstrom = 3E18 # Speed of light in angstroms
    
    if photometric_system == 'AB':
        AB_zero_point = 48.60 # 3161 Jy  after de-logging
        numer = np.sum( data * response * dlambda )
        denom = np.sum( (c_angstrom/(wavelength**2)) * response * dlambda )
        m_AB = -2.5*np.log10( numer/denom ) - AB_zero_point 
        return m_AB
    ##fi
#def

# ----------------------------------------------------------------------------

# First read the spectrum and the information file
spec_normed_dir = '../../data/spectra/eso_spectral_lib/normalized/'
spec_proper_dir = '../../data/spectra/eso_spectral_lib/proper/'

spec_library = ast1500.spectra.ESOSpectralLibrary(path_to_spectra=spec_normed_dir)
v_band_filter = ast1500.spectra.PhotometricFilter(filter_name='V',
    filter_class='Pickles')
ic_band_filter = ast1500.spectra.PhotometricFilter(filter_name='Ic',
    filter_class='Pickles')
k_band_filter = ast1500.spectra.PhotometricFilter(filter_name='K',
    filter_class='Pickles')

spec_info_filename = '../../data/spectra/eso_spectral_lib/information.fit'
hdu = fits.open(spec_info_filename)
data = hdu[1].data

# Get information about the spectrum
Mv_sptype = data['SpType']
Mv = data['Mbol']-data['BCv']
Mic = data['Mbol']-data['BCi']
Mk = data['Mbol']-data['BCk']

Mv_AB_correction = 0.02
Mic_AB_correction = 0.45
Mk_AB_correction = 1.85

# Loop over each spectrum
for i in range( spec_library.n_spectra ):
    
    # Read in the spectrum
    spec_name = spec_library.spectra_names[i][2:]
    spec_filename = spec_library.spectra_filenames[i]
    spec_wavelength, spec_data = np.genfromtxt(spec_filename,usecols=(0,1)).T
    
    # Figure out where this spectrum is in the information file
    where_spectra_in_info = np.where( Mv_sptype.upper() == spec_name.upper() )[0]
    this_Mv = Mv[where_spectra_in_info]
    this_Mic = Miv[where_spectra_in_info]
    this_Mk = Mk[where_spectra_in_info]
    
    # Calculate the filter responses and difference
    v_band_response = v_band_filter.response(spec_wavelength)
    ic_band_response = ic_band_filter.response(spec_wavelength)
    k_band_response = k_band_filter.response(spec_wavelength)
    delta_lambda = np.diff(spec_wavelength)[0]
    
    # Calculate the normalized V-band magnitude
    m_v_normed = calculate_magnitude(spec_data,v_band_response,spec_wavelength,
        delta_lambda)
    Mv_correction = this_Mv - mag_normed
    
    # Correct the data
    spec_data_corrected = spec_data * 10**(Mv_correction/(-2.5))
    Mv_new = calculate_magnitude(spec_data_corrected,v_band_response,spec_wavelength,
        delta_lambda)
    
    new_spec_filename = os.path.abspath(spec_proper_dir)+'/'+os.path.basename(spec_filename)
    np.savetxt(new_spec_filename, 
        np.array([spec_wavelength,spec_data_corrected]).T, delimiter='  ', 
        header='wavelength [Angstrom], F_lambda' )
###i
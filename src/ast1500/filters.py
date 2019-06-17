# ----------------------------------------------------------------------------
#
# TITLE -
# AUTHOR - James Lane
# PROJECT -
# CONTENTS:
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
'''
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb, glob


## Plotting
from matplotlib import pyplot as plt

## Astropy
# from astropy.io import fits
# from astropy.coordinates import SkyCoord
# from astropy import table
# from astropy import units as apu

# ----------------------------------------------------------------------------

## Example class
class ESOSpectralLibrary:
    '''
    ESOSpectralLibrary:

    Class to manipulate the spectral library of
    https://www.eso.org/sci/facilities/paranal/decommissioned/isaac/tools/lib.html

    Args:
        required_arg (data type) - Some description
        not_required_arg (data type) - Some description [the default]
    '''

    def __init__(self,
                 path_to_spectra='../../../data/eso_spectral_lib/spectra'
                ):
        # Set the path to the data
        self.path_to_spectra = os.path.abspath(path_to_spectra)

        # First get the names of all the filters
        spectra_fnames = glob.glob(path_to_spectra+'*.dat')
        if len(spectra_fnames)==0:
            raise ValueError('Could not find any .dat files in {}'.format(path_to_spectra))
        ##fi

        self.n_spectra = len(spectra_fnames)
        self.spectra_filenames = np.array([os.path.abspath(x) for x in spectra_fnames])
        self.spectra_names = np.array([os.path.basename(x)[:-4] for x in spectra_fnames])

        spectra_class = np.empty(self.n_spectra,dtype='U1')
        spectra_subclass = np.empty(self.n_spectra,dtype='U2')
        spectra_metallicity = np.empty(self.n_spectra,dtype='U1')
        spectra_luminosity = np.empty(self.n_spectra,dtype='U4')

        for i in range(self.n_spectra):
            temp_name = self.spectra_names[i].replace('uk','') # Remove uk

            # Check for metallicity information
            if temp_name[0] in ['w','r']:
                spectra_metallicity[i] = temp_name[0]
                temp_name = temp_name[1:]
            else:
                spectra_metallicity[i] = 'n'
            ##ie

            # Get spectral classification and luminosity type
            spectra_class[i]= temp_name[0]
            temp_name = temp_name[1:]
            spectra_subclass[i] = "".join([s for s in temp_name if s.isdigit()])
            spectra_luminosity[i] = "".join([s for s in temp_name if not s.isdigit()])
        ###i

        self.spectra_metalicity = spectra_metallicity
        self.spectra_class = spectra_class
        self.spectra_subclass = spectra_subclass
        self.spectra_luminosity = spectra_luminosity
    #def

    def spectra_class_available(self):
        '''spectra_classes_available

        Return the available spectral classes
        '''
        return np.unique(self.spectra_class)
    #def

    def find_spectra_of_class(self,spectra_class=None):
        '''find_spectra_of_class:

        Finds all the spectra of a given class and returns their names
        '''
        where_spectra_of_class = np.where( self.spectra_class == spectra_class )[0]
        return self.spectra_names[where_spectra_of_class]
    #def

    def _get_spectra_fname(self,name):
        if type(name) == int:
            fname = self.spectra_filenames[name]
        elif name in self.spectra_names:
            where_fname = np.where( self.spectra_names == name )[0]
            fname = self.spectra_filenames[where_fname]
        else:
            raise ValueError("name: {} cannot be found in the spectral library".format(name))
        ##ie
        return fname
    #def

    def read_spectra(self,name):
        '''read_spectra:

        Read a spectrum given an index number or the name of a spectrum. Return
        the raw data to the user.
        '''
        fname = self._get_spectra_fname(name)
        spectra_wavelength, spectra_data = np.genfromtxt(fname,usecols=(0,1)).T
        return spectra_wavelength, spectra_data
    #def

    def set_spectra(self,name):
        '''set_spectra:

        Set the current spectra
        '''
        fname = self._get_spectra_fname(name)
        spectra_wavelength, spectra_data = np.genfromtxt(fname,usecols=(0,1)).T
        self._cur_data_ = spectra_data
        self._cur_wavelength_ = spectra_wavelength
    #def

    def get_wavelength_range(self,name):
        '''get_wavelength_range:

        Find the wavelength range of the current spectra
        '''

#cls

# ----------------------------------------------------------------------------

class ESOSpectrum:
    '''ESOSpectrum

    Class to accesss ESO spectra
    '''

    def __init__(self,
                 fname=None,
                 data=None,
                 wavelength=None):
        # Two options to initialize: filename or data and wavelength
        pass
    #def
#cls

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

## Scipy
from scipy import interpolate

# ----------------------------------------------------------------------------

## Example class
class ESOSpectralLibrary:
    '''
    ESOSpectralLibrary:

    Class to manipulate the spectral library of
    https://www.eso.org/sci/facilities/paranal/decommissioned/isaac/tools/lib.html
    
    Data from Tables of the paper
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/PASP/110/863

    Args:
        required_arg (data type) - Some description
        not_required_arg (data type) - Some description [the default]
    '''

    def __init__(self,
                 path_to_spectra='/Users/JamesLane/Science/Projects/PhD/AST1500/data/spectra/eso_spectral_lib/proper'
                ):
        # Set the path to the data
        self.path_to_spectra = os.path.abspath(path_to_spectra)

        # First get the names of all the filters
        spectra_fnames = glob.glob(self.path_to_spectra+'/*.dat')
        if len(spectra_fnames)==0:
            raise ValueError('Could not find any .dat files in {}'.format(self.path_to_spectra))
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

        self.spectra_metallicity = spectra_metallicity
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
            where_fname = np.where( self.spectra_names == name )[0][0]
            fname = self.spectra_filenames[where_fname]
        else:
            raise ValueError("name: {} cannot be found in the spectral library".format(name))
        ##ie
        return fname
    #def

    def read_spectra(self,name,return_data=False,normalized=False):
        '''read_spectra:

        Read a spectrum given an index number or the name of a spectrum. Return
        the raw data to the user.
        '''
        fname = self._get_spectra_fname(name)
        if return_data:
            spectra_wavelength, spectra_data = np.genfromtxt(fname,usecols=(0,1)).T
            return spectra_wavelength, spectra_data
        else:
            spec = ESOSpectrum(filename=fname,normalized=normalized)
            return spec
        ##ie
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

    Class to accesss ESO spectra. Can be queried by either filename or 
    data and wavelength.
    
    Args:
        filename (str) -  
        data (float array) - 
        wavelength (float array) - 
        data_normalized (bool) - Has data input already been normalized 
            Such that it is f_lambda [False]
        normalization_method (int) - Method for normalizing spectral data. 
        
    '''

    def __init__(self,
                 filename=None,
                 data=None,
                 wavelength=None,
                 normalized=False
                 ):
                 
        # Two options to initialize: filename or data and wavelength
        self.normalized=normalized
        if filename == None:
            assert data != None and wavelength != None,\
                'If filename is not supplied then data and wavelength must be'
            self.filename = filename
            self.wavelength = wavelength
            if self.normalized:
                self.data = data
            else:
                self._data_raw = data
            ##ie
        elif filename != None:
            self.filename = filename
            read_wavelength,read_data_raw = np.genfromtxt(filename,usecols=(0,1)).T
            self.wavelength = read_wavelength
            if self.normalized:
                self.data = read_data_raw
            else:
                self._data_raw = read_data_raw
            ##ie
        ##fi  
        if self.normalized == False:
            print('Warning: spectra is not normalized')  
    ##fi
    
    def calculate_magnitude(self,photometric_filter=None,
        photometric_filter_kws={}):
        
        if photometric_filter != None:
            assert isinstance(photometric_filter,PhotometricFilter),\
                'photometric_filter must be a PhotometricFilter object'
        else:
            assert isinstance(photometric_filter_kws,dict),\
                'photometric_filter_kws must be a dictionary'
            assert len(photometric_filter_kws) > 0,\
                'photometric_filter_kws must have at least 1 element'
            
            photometric_filter = PhotometricFilter(**photometric_filter_kws)
        
        if self.normalized:
            response = photometric_filter.response(self.wavelength/10.)
            mag = -2.5*np.log10( self._convolve_filter_flux(response) )
        else:
            raise RuntimeError('Magnitudes for non-normalized spectra are not physical')
        ##ie
        return mag
    #def
    
    def plot_spectrum(self,fig=None,ax=None,range=None,normalize_to_1=True,plot_kws={}):
        '''plot_spectrum:
        
        Plot the spectrum of the source
        '''
        if range != None:
            assert len(range) == 2,'range must be a 2-element array'
            assert range[0] < range[1],\
                'The first element of range is not less than the second'
            where_plot_spectrum = np.where( ( self.wavelength > range[0]  ) &
                                            ( self.wavelength < range[1] ) )[0]
            plot_wavelength = self.wavelength(where_plot_spectrum)
            plot_data = self.data(where_plot_spectrum)
        ##fi
        if normalize_to_1 == True:
            plot_data = plot_data / np.max(plot_data)
        ##fi
        if fig == None or ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ##fi
        ax.plot(plot_wavelength,plot_data,**plot_kws)
        ax.set_xlabel(r'$\lambda [\AA]$')
        ax.set_ylabel(r'Flux')
        
        return fig,ax
    #def
    
    def _convolve_filter_flux(self,response):
        '''Assumes wavelength in angstroms
        '''
        c_angstrom = 3E18
        assert np.all(np.diff(self.wavelength)==np.diff(self.wavelength)[0]),\
            'Wavelength data must have equal spacing'
        dlambda=np.diff(self.wavelength)[0]
        return np.sum( self.data*response*dlambda ) /\
            np.sum( response*c_angstrom*dlambda/(self.wavelength**2) )
    ##def
#cls

# ------------------------------------------------------------------------------

class PhotometricFilter:
    '''PhotometricFilter

    Class to hold information about photometric filters. Currently supports:

    Bessell: U,B,V,R,I
    2MASS: J,H,Ks
    
    Can call just based on the name as long as there is no degenerate overlap. 
    By default the path is local to the module. If 
    
    Args:
        filter_name (string) - Name of the filter
        filter_name (string) - Class of the filter (e.g. Bessell, 2MASS, ..) [None]
        filter_path (string) - path to the filter information, None defaults
            to local module files [None]
        filter_response_wavelength (float array) - User-supplied filter response 
            wavelengths [None]
        filter_response_data (float array) - User-supplied filter response [None]
    '''
    def __init__(self,
                 filter_name,
                 filter_class=None,
                 filter_path=None,
                 filter_file_type='.dat',
                 filter_response_wavelength=None,
                 filter_response_data=None
                 ):
        # First determine if the filter will be read in or if the wavelength 
        # and response information will be provided as arguments
        if filter_response_data != None and filter_response_wavelength != None:
            _read_in_filter = False
        else:
            _read_in_filter = True
        ##ie
        
        # Either read in the filter or the use supplied the information
        if _read_in_filter:
        
            # Default path to filter data. It's a part of the module
            if filter_path == None:
                filter_path = os.path.dirname(__file__)+'/data/filters/'
            ##fi
            
            # As long as their is no overlap then set the class based just off the 
            # name of the filter
            if filter_name in ['U','V','B','R','I'] and filter_class == None:
                filter_class = 'Bessell'
            ##fi
            if filter_name in ['J','H','Ks'] and filter_class == None:
                filter_class = '2MASS'
            ##fi
            
            # Check to make sure the class and name make sense
            if filter_class not in ['Bessell','2MASS']:
                raise ValueError('{} filter class not supported'.format(filter_class))
            ##fi
            if filter_class == 'Bessell':
                if filter_name not in ['U','V','B','R','I']:
                    raise ValueError('{} not a Bessell filter'.format(filter_name))
                ##fi
            ##fi
            if filter_class == '2MASS':
                if filter_name not in ['J','H','Ks']:
                    raise ValueError('{} not a 2MASS filter'.format(filter_name))
                ##fi
            ##fi

            # Filter path stuff
            if filter_file_type[0] != '.':
                filter_file_type = '.'+filter_file_type
            ##fi
            filter_path = os.path.abspath(filter_path)
            filter_name_list = glob.glob(filter_path+'/*'+filter_file_type)
            if len(filter_name_list) == 0:
                raise RuntimeError('No filters of type '+filter_file_type+' in '+filter_path)
            ##fi
            assert filter_path+'/'+filter_class+'.'+filter_name+filter_file_type in filter_name_list,\
                filter_class+'.'+filter_name+filter_file_type+' not in '+filter_path
            ##as

            self.filter_name = filter_name
            self.filter_class = filter_class
            self.filename = filter_path+'/'+filter_class+'.'+filter_name+filter_file_type

            # Now read the filter
            filter_wavelength,filter_response = np.genfromtxt(self.filename).T
            self.wavelength_data = filter_wavelength/10. # In nm
            self.response_data = filter_response
        
        else:
            self.wavelength_data = filter_wavelength_data # In nm
            self.response_data = filter_response_data
            self.filter_class = filter_class # None unless supplied
            self.filter_name = filter_name # None unless supplied
        ##ie

        # Now construct a cubic spline of the filter response
        spl = interpolate.CubicSpline(self.wavelength_data,self.response_data)
        self._spline_ = spl
    #def

    def response(self,wavelength):
        '''response

        Use the cubic spline to get the filter response at an arbitrary
        wavelength

        Args:
            wavelength (float or array) - Wavelength(s) in nm to check response

        Returns:
            response (array) - Response of the filter at the wavelength
                (from 0 to 1)
        '''
        # Check if in wavelength range and output response
        wvln_min,wvln_max = self.get_wavelength_range()
        response = np.zeros_like(wavelength)
        where_in_filter_range = np.where( (wavelength > wvln_min) &
                                          (wavelength < wvln_max) )[0]
        response[where_in_filter_range] = self._spline_(wavelength[where_in_filter_range])
        return response
    #def

    def get_wavelength_range(self):
        '''get_wavelength_range:

        Get the wavelength range of the filter data

        Returns: wavelength_range (2-array) - wavelength minimum and
            maximum in nm
        '''
        return [np.min(self.wavelength_data),np.max(self.wavelength_data)]
    #def
#cls

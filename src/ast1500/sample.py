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
import sys, os, pdb
# import copy
# import glob
# import subprocess

## Plotting
from matplotlib import pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import colors
# from matplotlib import cm

## Astropy
# from astropy.io import fits
# from astropy.coordinates import SkyCoord
# from astropy import table
# from astropy import units as apu

## galpy
# from galpy import orbit
# from galpy import potential
# from galpy import actionAngle as aA
# from galpy import df
# from galpy.util import bovy_coords as gpcoords
# from galpy.util import bovy_conversion as gpconv
# from galpy.util import bovy_plot as gpplot

# ----------------------------------------------------------------------------

## Example class
class example_class:
    '''
    example:

    A description of the class

    Args:
        required_arg (data type) - Some description
        not_required_arg (data type) - Some description [the default]
    '''
    
    def __init__(   self,
                    required_arg,
                    not_required_arg=1.0
                ):
        self.property = required_arg
        self.other_property = not_required_arg
    #def
#cls

def example_function(   required_arg,
                        not_required_arg=1.0
                    ):
    '''
    example_function:
    
    Description of example function.
    
    Args:
        required_arg (data_type) - some description
        not_required_arg (data_type) - some description [default]
    
    Returns:
        some_variable (type) - some description
    
    Output:
        some_output (type) - some description
    
    Raises:
        TypeOfError - Some error that is caused by something
    '''
#def
    pass

# ----------------------------------------------------------------------------

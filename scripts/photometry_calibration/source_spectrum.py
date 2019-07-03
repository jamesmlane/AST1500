# ----------------------------------------------------------------------------
#
# TITLE -
# AUTHOR - James Lane
# PROJECT -
# CONTENTS:
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
''' Try and create a custom power law that does a good job of generating the 
observed zero point photon counts.
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os, pdb
from ipdb import set_trace

## Plotting
from matplotlib import pyplot as plt

## Astropy
from astropy import table

## Scipy
from scipy.optimize import fmin

# ----------------------------------------------------------------------------

## Read data
tab = table.Table.read('bandpass_defaults.txt', format='ascii.csv')

bp_names = tab['# BP name'].data
bp_cen = tab['BP center'].data
bp_width = tab['BP width'].data
bp_zp = tab['nPhotons'].data
n_bp = len(bp_names)

# Cut out certain bands
remove_bands = ['U','B','V0','V','R']
remove_inds = [ np.where(bp_names==band)[0][0] for band in remove_bands ]
bp_names = np.delete(bp_names, remove_inds)
bp_cen = np.delete(bp_cen, remove_inds)
bp_width = np.delete(bp_width, remove_inds)
bp_lambda_low = bp_cen - bp_width/2
bp_lambda_hi = bp_cen + bp_width/2
bp_zp = np.delete(bp_zp, remove_inds)
n_bp = len(bp_names)
data_arr = np.array([bp_lambda_low,bp_lambda_hi,bp_zp]).T
# data_arr = []
# for i in range(n_bp):
#     data_arr.append( [ bp_lambda_low[i], bp_lambda_hi[i], bp_zp[i] ] )
# ###i

# ----------------------------------------------------------------------------

# Now come up with a function that can be fit using Newtons method
def power_law(x,A,alpha,x0):
    return A*((x-x0)**alpha)
#def

def integrate_power_law(A,alpha,x0,lambda1,lambda2):
    '''integrate_power_law
    
    Integrate a power law between two wavelength cutoffs
    '''
    return (A/(alpha+1))*((lambda2-x0)**(alpha+1)-(lambda1-x0)**(alpha+1))
#def

def power_law_filter_sum_square(parms,x=data_arr):
    '''PL_filter
    
    x is an Nx3 vector of lower wavelength, upper wavelength, and integrated sum
    A is the amplitude of the power law profile
    alpha is the power law index
    
    Calculates the squared sum of the difference between the expected
    integrated number of photons and the calculated integrated number of
    photons.
    '''
    A,alpha,x0 = parms
    
    sum_abs_difference = 0
    for i in range(len(x)):
        lambda1,lambda2,known_nphotons = x[i]
        # Integrate the power law between the two limits
        calc_nphotons = integrate_power_law(A,alpha,x0,lambda1,lambda2)
        sum_abs_difference += np.square( calc_nphotons - known_nphotons )
    ###i    
    return sum_abs_difference
#def

def plot_power_law(data_arr, A, alpha, x0, names):
    '''plot_power_law:
    
    '''
    fig = plt.figure(figsize=(8,8))
    axs = fig.subplots(ncols=1, nrows=2)
    
    bp_lambda_low, bp_lambda_hi, bp_zp = data_arr.T
    
    
    # First axis
    mock_x = np.arange( np.amin(bp_lambda_low)-1e-7, np.amax(bp_lambda_hi)+1e-7, 1e-8 )
    axs[0].plot( mock_x*1e6, power_law(mock_x,A,alpha,x0), color='Black' )
    axs[0].set_xlabel(r'$\lambda$ [$\mu$m]')
    axs[0].set_ylabel(r'$N_{gamma}$')
    
    # Second axis
    assert len(bp_names) == len(bp_zp)
    for i in range(len(bp_names)):
        bp_cen = bp_lambda_low[i] + ( bp_lambda_hi[i] - bp_lambda_low[i] )/2
        calc_np = integrate_power_law(A,alpha,x0,bp_lambda_low[i],bp_lambda_hi[i])
        perc_diff = 100*((calc_np-bp_zp[i])/bp_zp[i])
        axs[1].scatter( bp_cen*1e6, perc_diff, color='Black' )
        axs[1].annotate( bp_names[i], xy=(bp_cen*1e6, perc_diff*1.05),
        xytext=(0.1+bp_cen*1e6, perc_diff) )
    ###i
    
    axs[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    axs[1].set_ylabel(r'$N_{\gamma,calc}-N_{\gamma,BP}$ [% diff.]')
    
    return fig, axs
#def

# ----------------------------------------------------------------------------

# result = fmin( power_law_filter_sum_square, x0=[1,-3.5,1e-7] )
# print(result)
result = [ 1.23341197e-01, -3.45115552e+00, -6.49788953e-07]

fig, axs = plot_power_law( data_arr, result[0], result[1], result[2], bp_names )

fig.savefig('test.png')

# set_trace()        
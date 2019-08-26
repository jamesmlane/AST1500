# ----------------------------------------------------------------------------
#
# TITLE - photometry_calibration.py
# AUTHOR - James Lane
# PROJECT - AST 1500
#
# ----------------------------------------------------------------------------

### Docstrings and metadata:
''' Investigate the zero-points for OOMAO's photometry class
'''
__author__ = "James Lane"

### Imports

## Basic
import numpy as np
import sys, os
# from ipdb import set_trace
import pdb

## Plotting
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import colors
# from matplotlib import cm

## Astropy
from astropy import table

# Scipy
from scipy.optimize import curve_fit

# ----------------------------------------------------------------------------

## Read data
tab = table.Table.read('bandpass_defaults.txt', format='ascii.csv')

bp_names = tab['# BP name'].data
bp_cen = tab['BP center'].data
bp_width = tab['BP width'].data
bp_zp = tab['nPhotons'].data
n_bp = len(bp_names)

# Cut out certain bands
remove_bands = ['U','B','V0','V']
remove_inds = [ np.where(bp_names==band)[0][0] for band in remove_bands ]
bp_names = np.delete(bp_names, remove_inds)
bp_cen = np.delete(bp_cen, remove_inds)
bp_width = np.delete(bp_width, remove_inds)
bp_zp = np.delete(bp_zp, remove_inds)
n_bp = len(bp_names)

# ----------------------------------------------------------------------------

# fig = plt.figure()
# ax = fig.add_subplot(111)
# 
# for i in range(n_bp):
#     ax.scatter( bp_cen[i]*1e6, bp_zp[i]/1e11, label=bp_names[i] )
#     ax.plot([(bp_cen[i]-bp_width[i])*1e6,(bp_cen[i]+bp_width[i])*1e6],
#             [bp_zp[i]/1e11,bp_zp[i]/1e11])
# 
# ax.set_xlabel('Bandpass ($\mu$m)')
# ax.set_ylabel('Zero point ($10^{11}$ photons)')
# ax.legend()
# 
# plt.savefig('bandpass_plots.pdf')

# ----------------------------------------------------------------------------

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# 
# for i in range(n_bp):
#     ax.scatter( bp_cen[i]*1e6, bp_width[i]*1e6, bp_zp[i]/1e11 )
#     ax.text( bp_cen[i]*1e6, bp_width[i]*1e6, bp_zp[i]/1e11, 
#             '%s' % (bp_names[i]), size=20, zorder=1, color='Black')
# ###i
# 
# ax.set_xlabel('BP Center ($\mu$m)')
# ax.set_ylabel('BP Width ($\mu$m)')
# ax.set_zlabel('Zero point ($10^{11}$ photons)')
# 
# plt.show()

# ------------------------------------------------------------------------------

# Fit a 2-D power law surface to the zero-points

def power_law_surf_fit(x,A,B,k1,k2,d):
    return A*x[0]**k1 + B*x[1]**k2 + d
#def
def power_law_surf(x,y,A,B,k1,k2,d):
    return A*x**k1 + B*y**k2 + d
#def

fit_data = np.vstack((bp_cen,bp_width))

popt,pcov = curve_fit( power_law_surf_fit, fit_data*1e6, bp_zp/1e11, maxfev=100000 )

fit_zp = power_law_surf_fit(fit_data*1e6,*popt)

z_wc = 9500 / 1E4 # In microns
Y_wc = 10200  / 1E4 # In microns
z_bp = (10700-8200) / 1E4 # In microns
Y_bp = (11000-9400) / 1E4 # In microns

z_zp = power_law_surf_fit([z_wc,z_bp],*popt)
Y_zp = power_law_surf_fit([Y_wc,Y_bp],*popt)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(n_bp):
    ax.scatter( bp_cen[i]*1e6, bp_width[i]*1e6, bp_zp[i]/1e11 )
    ax.scatter( bp_cen[i]*1e6, bp_width[i]*1e6, fit_zp[i], color='Black' )
    ax.text( bp_cen[i]*1e6, bp_width[i]*1e6, bp_zp[i]/1e11, 
            '%s' % (bp_names[i]), size=20, zorder=1, color='Black')
###i

ax.scatter( z_wc, z_bp, z_zp, color='Black' )
ax.text( z_wc, z_bp, z_zp, 
        '%s' % ('z'), size=20, zorder=1, color='Black')
ax.scatter( Y_wc, Y_bp, Y_zp, color='Black' )
ax.text( Y_wc, Y_bp, Y_zp, 
        '%s' % ('Y'), size=20, zorder=1, color='Black')

print('z ZP: '+str(z_zp))
print('Y ZP: '+str(Y_zp))

ax.set_xlabel('BP Center ($\mu$m)')
ax.set_ylabel('BP Width ($\mu$m)')
ax.set_zlabel('Zero point ($10^{11}$ photons)')

plt.show()
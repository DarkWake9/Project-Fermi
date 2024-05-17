import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.io.fits as pyfits
from astropy.table import Table
from astropy import units as u
from scipy import interpolate
from scipy.integrate import quad
dir = '3B-20180206-20180622/'

lc_file_path = 'Output/Light_curve_001/4fgl_j0854.8+2006_lightcurve.fits'
emin = 0.1 *u.GeV
emax = 300 *u.GeV
D_l = 134.1 *u.Mpc
hdul = fits.open('./' + dir+lc_file_path)
lc = hdul[1].data
Table(lc).colnames
lc['dnde']
#### Plotting spline fit
# flux = lc[1]['flux']
# flux_err = lc[1]['flux_err']

# tmean = (lc['tmin_mjd'] + lc['tmax_mjd'])/2
# TSmin = 9

# if len(lc['eflux'][lc['ts']>TSmin]) > 0:
#     scale = int(np.log10(lc['eflux'][lc['ts']>TSmin].max())) -2
# else:
#     scale = int(np.log10(lc['eflux_ul95'].max())) -2

# f = plt.figure(figsize=(10, 5))
# ax = f.add_subplot(111)
# ax.tick_params(which='major', length=5, direction='in')
# ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
# ax.tick_params(bottom=True, top=True, left=True, right=True)

# # if len(tmean[lc['ts']>TSmin]) > 9:
# time_continuum = np.linspace(np.min(lc['tmin_mjd']),np.max(lc['tmax_mjd']),10*len(lc['tmin_mjd']))
# tck_flux = interpolate.splrep(tmean[lc['ts']>TSmin], (10**-scale)*lc['flux'][lc['ts']>TSmin], k=3)
# tck_flux_error = interpolate.splrep(tmean[lc['ts']>TSmin], (10**-scale)*lc['flux_err'][lc['ts']>TSmin],k=3)

# flux_continuum = interpolate.splev(time_continuum, tck_flux)
# flux_continuum_err = interpolate.splev(time_continuum, tck_flux_error)
# ax.plot(time_continuum,flux_continuum,color="C1", label="Spline")



# plt.errorbar(tmean[lc['ts']>TSmin], (10**-scale)*lc['flux'][lc['ts']>TSmin], xerr = [ tmean[lc['ts']>TSmin]- lc['tmin_mjd'][lc['ts']>TSmin], lc['tmax_mjd'][lc['ts']>TSmin] - tmean[lc['ts']>TSmin] ], yerr=(10**-scale)*lc['flux_err'][lc['ts']>TSmin], markeredgecolor='black', fmt='o', capsize=4)
# plt.errorbar(tmean[lc['ts']<=TSmin], (10**-scale)*lc['flux_ul95'][lc['ts']<=TSmin], xerr = [ tmean[lc['ts']<=TSmin]- lc['tmin_mjd'][lc['ts']<=TSmin], lc['tmax_mjd'][lc['ts']<=TSmin] - tmean[lc['ts']<=TSmin] ], yerr=5*np.ones(len(lc['flux_err'][lc['ts']<=TSmin])), markeredgecolor='black', fmt='o', uplims=True, color='orange', capsize=4)
# plt.xlabel('Time [MJD]')
# plt.ylabel('Flux * 1e-9 [$erg$ $cm^{-2}$ $s^{-1}$]')

# if len(lc['flux'][lc['ts']>TSmin]) > 0:
#     y0 = (lc['flux'][lc['ts']>TSmin]).max()
#     y1 = (lc['flux'][lc['ts']>TSmin] + lc['flux_err'][lc['ts']>TSmin]).max()
#     if y1 > 4*y0:
#         y1 = 4*y0
    
#     if len(lc['flux_ul95'][lc['ts']<=TSmin]) > 0:
#         y2 = (lc['flux_ul95'][lc['ts']<=TSmin]).max()
#         if y2 > y1:
#             y1 = y2
    
# else:
#     y1 = (lc['flux_ul95'][lc['ts']<=TSmin]).max()



# ymin = -(10**-scale)*0.1*y1               
# plt.ylim(ymin,(10**-scale)*1.1*y1)
# plt.legend()
# plt.title('Light curve')

# ax.grid(linestyle=':',which='both')
# plt.show()
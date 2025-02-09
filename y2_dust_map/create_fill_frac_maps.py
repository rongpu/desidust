from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/imaging_systematics'))
from plot_healpix_map import plot_map


v_missing = 1.
v_good = 0.

fwhm_dict = {128: 0.5, 256: 0.25, 512: 0.125}  # degree
n_iteration_dict = {128: 200, 256: 400, 512: 800}

# nside = 128
# fwhm = 0.5
# n_iteration = 200

# nside = 256
# fwhm = 0.25  # degree
# n_iteration = 400

# nside = 512
# fwhm = 0.125  # degree
# n_iteration = 800


################################################################################################################################################

for nside in [128, 256, 512]:

    fwhm = fwhm_dict[nside]
    n_iteration = n_iteration_dict[nside]

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/dgr_map_combined_{}.fits'.format(nside)))
    maps['FILL_FRAC'] = v_good
    maps = maps[['HPXPIXEL', 'FILL_FRAC']]

    tmp = maps.copy()
    tmp.sort('HPXPIXEL')
    pixels, v = tmp['HPXPIXEL'], tmp['FILL_FRAC']

    npix = hp.nside2npix(nside)
    map_values = np.full(npix, v_missing)
    map_values[pixels] = v
    hp_mask = np.in1d(np.arange(npix), pixels)
    newmap = hp.ma(map_values)
    newmap.data[~hp_mask] = v_missing

    for ii in range(n_iteration):
        smoothed_map = hp.sphtfunc.smoothing(newmap, fwhm=fwhm/(180/np.pi))
        newmap.data[~hp_mask] = smoothed_map[~hp_mask]

    maps_fill_frac = Table()
    maps_fill_frac['HPXPIXEL'] = np.arange(len(smoothed_map))
    maps_fill_frac['FILL_FRAC'] = smoothed_map.data

    maps_fill_frac.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1_fill_frac_gr_{}.fits'.format(nside), overwrite=True)

    default_dpi = {32: 100, 64: 200, 128: 400, 256: 600, 512: 1200}
    default_xsize = {32: 1500, 64: 4000, 128: 4000, 256: 6000, 512: 12000}
    plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/dust/desi_ebv/v1/tmp'

    fn = os.path.join(plot_dir, 'v1_fill_frac_{}_gr.png'.format(nside))
    plot_map(nside, maps_fill_frac['FILL_FRAC'], maps_fill_frac['HPXPIXEL'],
             vmin=0., vmax=1., cmap='jet', dpi=default_dpi[nside], xsize=default_xsize[nside],
             cbar_label='$E(B-V)_\mathrm{DESI}$ (mag)', save_path=fn, show=False, timing=False)
    # display(Image(fn))

################################################################################################################################################

for nside in [128, 256, 512]:

    fwhm = fwhm_dict[nside]
    n_iteration = n_iteration_dict[nside]

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_combined_{}.fits'.format(nside)))
    maps['FILL_FRAC'] = v_good
    maps = maps[['HPXPIXEL', 'FILL_FRAC']]

    tmp = maps.copy()
    tmp.sort('HPXPIXEL')
    pixels, v = tmp['HPXPIXEL'], tmp['FILL_FRAC']

    npix = hp.nside2npix(nside)
    map_values = np.full(npix, v_missing)
    map_values[pixels] = v
    hp_mask = np.in1d(np.arange(npix), pixels)
    newmap = hp.ma(map_values)
    newmap.data[~hp_mask] = v_missing

    for ii in range(n_iteration):
        smoothed_map = hp.sphtfunc.smoothing(newmap, fwhm=fwhm/(180/np.pi))
        newmap.data[~hp_mask] = smoothed_map[~hp_mask]

    maps_fill_frac = Table()
    maps_fill_frac['HPXPIXEL'] = np.arange(len(smoothed_map))
    maps_fill_frac['FILL_FRAC'] = smoothed_map.data

    maps_fill_frac.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1_fill_frac_rz_{}.fits'.format(nside), overwrite=True)

    default_dpi = {32: 100, 64: 200, 128: 400, 256: 600, 512: 1200}
    default_xsize = {32: 1500, 64: 4000, 128: 4000, 256: 6000, 512: 12000}
    plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/dust/desi_ebv/v1/tmp'

    fn = os.path.join(plot_dir, 'v1_fill_frac_{}_rz.png'.format(nside))
    plot_map(nside, maps_fill_frac['FILL_FRAC'], maps_fill_frac['HPXPIXEL'],
             vmin=0., vmax=1., cmap='jet', dpi=default_dpi[nside], xsize=default_xsize[nside],
             cbar_label='$E(B-V)_\mathrm{DESI}$ (mag)', save_path=fn, show=False, timing=False)
    # display(Image(fn))

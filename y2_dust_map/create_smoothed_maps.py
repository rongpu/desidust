from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/imaging_systematics'))
from plot_healpix_map import plot_map


fwhm_dict = {128: 0.5, 256: 0.25, 512: 0.125}  # degree
n_iteration_dict = {128: 200, 256: 400, 512: 800}

################################################################################################################################################

for nside in [128, 256, 512]:

    fwhm = fwhm_dict[nside]
    n_iteration = n_iteration_dict[nside]

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/dgr_map_combined_{}.fits'.format(nside)))
    maps['delta_gr'] = maps['delta_gr_wmean'].copy()
    maps = maps[['HPXPIXEL', 'delta_gr', 'n_star', 'EBV_SFD']]

    tmp = maps.copy()
    tmp.sort('HPXPIXEL')
    pixels, v = tmp['HPXPIXEL'], tmp['delta_gr']

    v_fill = np.median(v)

    npix = hp.nside2npix(nside)
    map_values = np.full(npix, v_fill)
    map_values[pixels] = v
    hp_mask = np.in1d(np.arange(npix), pixels)
    newmap = hp.ma(map_values)
    newmap.data[~hp_mask] = v_fill

    # plot_map(nside, newmap,
    #          cmap='jet', vmin=v_fill-0.02, vmax=v_fill+0.02,
    #          title=' ', show=True, xsize=1500, dpi=100)

    for ii in range(n_iteration):
        smoothed_map = hp.sphtfunc.smoothing(newmap, fwhm=fwhm/(180/np.pi))
        newmap.data[~hp_mask] = smoothed_map[~hp_mask]
        # if ii%20==0:
        #     plot_map(nside, newmap,
        #              cmap='jet', vmin=v_fill-0.02, vmax=v_fill+0.02,
        #              title=' ', show=True, xsize=1500, dpi=100)

    maps_filled = Table()
    maps_filled['HPXPIXEL'] = np.arange(len(newmap))
    maps_filled['delta_gr'] = newmap.data
    maps_filled = join(maps_filled, maps[['HPXPIXEL', 'n_star']], keys='HPXPIXEL', join_type='left').filled(0)
    maps_filled.sort('HPXPIXEL')

    maps_smooth = Table()
    maps_smooth['HPXPIXEL'] = np.arange(len(smoothed_map))
    maps_smooth['delta_gr'] = smoothed_map.data

    maps_filled.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1.1_desi_delta_gr_filled_{}.fits'.format(nside), overwrite=True)
    maps_smooth.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1.1_desi_delta_gr_smooth_{}.fits'.format(nside), overwrite=True)

    default_dpi = {32: 100, 64: 200, 128: 400, 256: 600, 512: 1200}
    default_xsize = {32: 1500, 64: 4000, 128: 4000, 256: 6000, 512: 12000}
    plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/dust/desi_ebv/v1.1/tmp'

    fn = os.path.join(plot_dir, 'v1.1_ebv_desi_{}_gr_smooth.png'.format(nside))
    plot_map(nside, maps_smooth['delta_gr'], maps_smooth['HPXPIXEL'],
             vmin=0., vmax=0.2, cmap='jet', dpi=default_dpi[nside], xsize=default_xsize[nside],
             cbar_label='$E(B-V)_\mathrm{DESI}$ (mag)', save_path=fn, show=False, timing=False)
    # display(Image(fn))

    fn = os.path.join(plot_dir, 'v1.1_ebv_desi_{}_gr_filled.png'.format(nside))
    plot_map(nside, maps_filled['delta_gr'], maps_filled['HPXPIXEL'],
             vmin=0., vmax=0.2, cmap='jet', dpi=default_dpi[nside], xsize=default_xsize[nside],
             cbar_label='$E(B-V)_\mathrm{DESI}$ (mag)', save_path=fn, show=False, timing=False)
    # display(Image(fn))

################################################################################################################################################

for nside in [128, 256, 512]:

    fwhm = fwhm_dict[nside]
    n_iteration = n_iteration_dict[nside]

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_combined_{}.fits'.format(nside)))
    maps['delta_rz'] = maps['delta_rz_wmean'].copy()
    maps = maps[['HPXPIXEL', 'delta_rz', 'n_star', 'EBV_SFD']]

    tmp = maps.copy()
    tmp.sort('HPXPIXEL')
    pixels, v = tmp['HPXPIXEL'], tmp['delta_rz']

    v_fill = np.median(v)

    npix = hp.nside2npix(nside)
    map_values = np.full(npix, v_fill)
    map_values[pixels] = v
    hp_mask = np.in1d(np.arange(npix), pixels)
    newmap = hp.ma(map_values)
    newmap.data[~hp_mask] = v_fill

    # plot_map(nside, newmap,
    #          cmap='jet', vmin=v_fill-0.02, vmax=v_fill+0.02,
    #          title=' ', show=True, xsize=1500, dpi=100)

    for ii in range(n_iteration):
        smoothed_map = hp.sphtfunc.smoothing(newmap, fwhm=fwhm/(180/np.pi))
        newmap.data[~hp_mask] = smoothed_map[~hp_mask]
        # if ii%20==0:
        #     plot_map(nside, newmap,
        #              cmap='jet', vmin=v_fill-0.02, vmax=v_fill+0.02,
        #              title=' ', show=True, xsize=1500, dpi=100)

    maps_filled = Table()
    maps_filled['HPXPIXEL'] = np.arange(len(newmap))
    maps_filled['delta_rz'] = newmap.data
    maps_filled = join(maps_filled, maps[['HPXPIXEL', 'n_star']], keys='HPXPIXEL', join_type='left').filled(0)
    maps_filled.sort('HPXPIXEL')

    maps_smooth = Table()
    maps_smooth['HPXPIXEL'] = np.arange(len(smoothed_map))
    maps_smooth['delta_rz'] = smoothed_map.data

    maps_filled.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1.1_desi_delta_rz_filled_{}.fits'.format(nside), overwrite=True)
    maps_smooth.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1.1_desi_delta_rz_smooth_{}.fits'.format(nside), overwrite=True)

    default_dpi = {32: 100, 64: 200, 128: 400, 256: 600, 512: 1200}
    default_xsize = {32: 1500, 64: 4000, 128: 4000, 256: 6000, 512: 12000}
    plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/dust/desi_ebv/v1.1/tmp'

    fn = os.path.join(plot_dir, 'v1.1_ebv_desi_{}_rz_smooth.png'.format(nside))
    plot_map(nside, maps_smooth['delta_rz'], maps_smooth['HPXPIXEL'],
             vmin=0., vmax=0.2, cmap='jet', dpi=default_dpi[nside], xsize=default_xsize[nside],
             cbar_label='$E(B-V)_\mathrm{DESI}$ (mag)', save_path=fn, show=False, timing=False)
    # display(Image(fn))

    fn = os.path.join(plot_dir, 'v1.1_ebv_desi_{}_rz_filled.png'.format(nside))
    plot_map(nside, maps_filled['delta_rz'], maps_filled['HPXPIXEL'],
             vmin=0., vmax=0.2, cmap='jet', dpi=default_dpi[nside], xsize=default_xsize[nside],
             cbar_label='$E(B-V)_\mathrm{DESI}$ (mag)', save_path=fn, show=False, timing=False)
    # display(Image(fn))



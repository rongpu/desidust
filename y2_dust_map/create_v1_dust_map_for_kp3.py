from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

import healpy as hp


# Create smoothed maps for filling in the missing pixels

for nside in [64, 128, 256, 512]:

    if nside==256 or nside==512:
        nside_smooth = nside
    else:
        nside_smooth = 128

    ##############################################################################################################

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/dgr_map_combined_{}.fits'.format(nside)))
    maps['EBV_DESI'] = maps['delta_gr_wmean']/1.049
    maps['EBV_DESI_ERR'] = maps['delta_gr_err']/1.049
    maps = maps[['HPXPIXEL', 'EBV_DESI', 'EBV_DESI_ERR', 'n_star', 'EBV_SFD']]
    maps_smooth = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1_desi_ebv_gr_smooth_{}.fits'.format(nside_smooth)))
    ebv = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_ring.npy'.format(nside_smooth))
    maps_smooth['EBV_SFD'] = ebv.copy()

    if nside!=nside_smooth:
        tmp = Table()
        tmp['HPXPIXEL'] = np.arange(hp.nside2npix(nside))
        tmp['EBV_DESI'] = hp.ud_grade(maps_smooth['EBV_DESI'], nside, order_in='RING', order_out='RING').data
        tmp['EBV_SFD'] = hp.ud_grade(maps_smooth['EBV_SFD'], nside, order_in='RING', order_out='RING').data
        maps_smooth = tmp.copy()
    maps_smooth['EBV_DESI_ERR'] = np.nan

    mask = ~np.in1d(maps_smooth['HPXPIXEL'], maps['HPXPIXEL'])
    maps_smooth = maps_smooth[mask]
    maps = vstack([maps, maps_smooth]).filled(0)
    maps.sort('HPXPIXEL')
    grmaps = maps.copy()

    ##############################################################################################################

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_combined_{}.fits'.format(nside)))
    maps['EBV_DESI'] = maps['delta_rz_wmean']/0.954
    maps['EBV_DESI_ERR'] = maps['delta_rz_err']/0.954
    maps = maps[['HPXPIXEL', 'EBV_DESI', 'EBV_DESI_ERR', 'n_star', 'EBV_SFD']]
    maps_smooth = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1_desi_ebv_rz_smooth_{}.fits'.format(nside_smooth)))
    maps_smooth['EBV_SFD'] = ebv.copy()

    if nside!=nside_smooth:
        tmp = Table()
        tmp['HPXPIXEL'] = np.arange(hp.nside2npix(nside))
        tmp['EBV_DESI'] = hp.ud_grade(maps_smooth['EBV_DESI'], nside, order_in='RING', order_out='RING').data
        tmp['EBV_SFD'] = hp.ud_grade(maps_smooth['EBV_SFD'], nside, order_in='RING', order_out='RING').data
        maps_smooth = tmp.copy()
    maps_smooth['EBV_DESI_ERR'] = np.nan

    mask = ~np.in1d(maps_smooth['HPXPIXEL'], maps['HPXPIXEL'])
    maps_smooth = maps_smooth[mask]
    maps = vstack([maps, maps_smooth]).filled(0)
    maps.sort('HPXPIXEL')
    rzmaps = maps.copy()

    ##############################################################################################################

    grmaps.rename_columns(['HPXPIXEL', 'EBV_DESI', 'EBV_DESI_ERR', 'n_star', 'EBV_SFD'], ['HPXPIXEL', 'EBV_DESI_GR', 'EBV_DESI_ERR_GR', 'N_STAR_GR', 'EBV_SFD_GR'])
    rzmaps.rename_columns(['HPXPIXEL', 'EBV_DESI', 'EBV_DESI_ERR', 'n_star', 'EBV_SFD'], ['HPXPIXEL', 'EBV_DESI_RZ', 'EBV_DESI_ERR_RZ', 'N_STAR_RZ', 'EBV_SFD_RZ'])

    maps = join(grmaps, rzmaps, keys='HPXPIXEL', join_type='outer')
    assert np.all(maps['HPXPIXEL']==np.arange(len(maps)))

    maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/kp3_maps/v1_desi_ebv_{}.fits'.format(nside), overwrite=False)

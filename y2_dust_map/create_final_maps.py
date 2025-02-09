# Create the final maps for public release

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from get_ebv_from_map import get_ebv_from_map


r = {'g_south': 3.214, 'r_south': 2.165, 'z_south': 1.211, 'g_north': 3.258, 'r_north': 2.176, 'z_north': 1.199}

fillfrac_thresholds = {128: 0.01, 256: 0.005, 512: 0.001}

################################# E(g-r) #################################
print('gr')

for nside in [128, 256, 512]:

    print('NSIDE', nside)

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/dgr_map_combined_{}.fits'.format(nside)))
    maps = maps[['HPXPIXEL', 'delta_gr_wmean', 'delta_gr_err', 'n_star', 'EBV_SFD']]
    maps.rename_columns(['delta_gr_wmean', 'delta_gr_err'], ['EGR', 'EGR_ERR'])

    maps_smooth = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1.1_desi_delta_gr_smooth_{}.fits'.format(nside)))
    ebv = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_ring.npy'.format(nside))
    maps_smooth['EBV_SFD'] = ebv.copy()
    maps_smooth.rename_columns(['delta_gr'], ['EGR'])

    mask = ~np.in1d(maps_smooth['HPXPIXEL'], maps['HPXPIXEL'])
    maps_smooth = maps_smooth[mask]
    maps = vstack([maps, maps_smooth]).filled(0)
    maps.sort('HPXPIXEL')

    maps['EBV_GR'] = maps['EGR'] / (r['g_south']-r['r_south'])
    maps['EBV_GR_ERR'] = maps['EGR_ERR'] / (r['g_south']-r['r_south'])

    nan_cols = ['EGR_ERR', 'EBV_GR_ERR']
    mask = maps['n_star']==0
    for col in nan_cols:
        maps[col][mask] = np.nan

    fillfrac = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/v1/smoothed/v1_fill_frac_gr_{}.fits'.format(nside)))
    assert np.all(maps['HPXPIXEL']==fillfrac['HPXPIXEL'])
    maps['FILL_FRAC'] = fillfrac['FILL_FRAC']

    mask_sky = maps['FILL_FRAC'] < fillfrac_thresholds[nside]
    mask_sky |= maps['n_star']>0
    fsky = np.sum(mask_sky)/len(mask_sky)
    print('Fraction of sky covered: {:.1f}%'.format(100*fsky))
    print('Fraction of empty pixels in the "good" area: {:.1f}%'.format(100*np.sum(mask_sky & (maps['n_star']==0))/np.sum(mask_sky)))

    for col in maps.colnames:
        maps.rename_column(col, col.upper())

    # Write "all-sky" map (not for public use)
    maps1 = maps.copy()
    maps1['GOOD'] = mask_sky.copy()
    maps1.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/desi_dust_gr_{}_allsky.fits'.format(nside), overwrite=False)

    maps = maps[mask_sky]

    maps['RA'], maps['DEC'] = hp.pix2ang(nside, maps['HPXPIXEL'], nest=False, lonlat=True)

    maps = maps[['HPXPIXEL', 'RA', 'DEC', 'N_STAR', 'EGR', 'EGR_ERR', 'EBV_GR', 'EBV_GR_ERR', 'EBV_SFD', 'FILL_FRAC']]

    # use int32 and float32 to reduce file size
    for col in maps.colnames:
        if col in ['HPXPIXEL', 'RA', 'DEC']:
            continue
        if maps[col].dtype=='>i8' or maps[col].dtype=='int64':
            maps[col] = np.array(maps[col], dtype='int32')
        elif maps[col].dtype=='>f8' or maps[col].dtype=='float64':
            maps[col] = np.array(maps[col], dtype='float32')

    output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/final_maps/desi_dust_gr_{}.fits'.format(nside)
    # maps.write(output_fn, overwrite=False)
    # Add header and extension names
    hdus = fits.HDUList()
    hdr = fits.Header()
    hdr['NSIDE'] = nside
    hdr['NEST'] = False
    hdus.append(fits.PrimaryHDU(data=None, header=hdr))
    hdus.append(fits.BinTableHDU(maps, name='MAPS'))
    hdus.writeto(output_fn, overwrite=False)

################################# E(r-z) #################################
print('rz')

for nside in [128, 256, 512]:

    print('NSIDE', nside)

    maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_combined_{}.fits'.format(nside)))
    maps = maps[['HPXPIXEL', 'delta_rz_wmean', 'delta_rz_err', 'n_star', 'EBV_SFD']]
    maps.rename_columns(['delta_rz_wmean', 'delta_rz_err'], ['ERZ', 'ERZ_ERR'])

    maps_smooth = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/v1.1_desi_delta_rz_smooth_{}.fits'.format(nside)))
    ebv = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_ring.npy'.format(nside))
    maps_smooth['EBV_SFD'] = ebv.copy()
    maps_smooth.rename_columns(['delta_rz'], ['ERZ'])

    mask = ~np.in1d(maps_smooth['HPXPIXEL'], maps['HPXPIXEL'])
    maps_smooth = maps_smooth[mask]
    maps = vstack([maps, maps_smooth]).filled(0)
    maps.sort('HPXPIXEL')

    maps['EBV_RZ'] = maps['ERZ'] / (r['r_south']-r['z_south'])
    maps['EBV_RZ_ERR'] = maps['ERZ_ERR'] / (r['r_south']-r['z_south'])

    nan_cols = ['ERZ_ERR', 'EBV_RZ_ERR']
    mask = maps['n_star']==0
    for col in nan_cols:
        maps[col][mask] = np.nan

    fillfrac = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/v1/smoothed/v1_fill_frac_rz_{}.fits'.format(nside)))
    assert np.all(maps['HPXPIXEL']==fillfrac['HPXPIXEL'])
    maps['FILL_FRAC'] = fillfrac['FILL_FRAC']

    mask_sky = maps['FILL_FRAC'] < fillfrac_thresholds[nside]
    mask_sky |= maps['n_star']>0
    fsky = np.sum(mask_sky)/len(mask_sky)
    print('Fraction of sky covered: {:.1f}%'.format(100*fsky))
    print('Fraction of empty pixels in the "good" area: {:.1f}%'.format(100*np.sum(mask_sky & (maps['n_star']==0))/np.sum(mask_sky)))

    for col in maps.colnames:
        maps.rename_column(col, col.upper())

    # Write "all-sky" map (not for public use)
    maps1 = maps.copy()
    maps1['GOOD'] = mask_sky.copy()
    maps1.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/desi_dust_rz_{}_allsky.fits'.format(nside), overwrite=False)

    maps = maps[mask_sky]

    maps['RA'], maps['DEC'] = hp.pix2ang(nside, maps['HPXPIXEL'], nest=False, lonlat=True)

    maps = maps[['HPXPIXEL', 'RA', 'DEC', 'N_STAR', 'ERZ', 'ERZ_ERR', 'EBV_RZ', 'EBV_RZ_ERR', 'EBV_SFD', 'FILL_FRAC']]

    # use int32 and float32 to reduce file size
    for col in maps.colnames:
        if col in ['HPXPIXEL', 'RA', 'DEC']:
            continue
        if maps[col].dtype=='>i8' or maps[col].dtype=='int64':
            maps[col] = np.array(maps[col], dtype='int32')
        elif maps[col].dtype=='>f8' or maps[col].dtype=='float64':
            maps[col] = np.array(maps[col], dtype='float32')

    output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/final_maps/desi_dust_rz_{}.fits'.format(nside)
    # maps.write(output_fn, overwrite=False)
    # Add header and extension names
    hdus = fits.HDUList()
    hdr = fits.Header()
    hdr['NSIDE'] = nside
    hdr['NEST'] = False
    hdus.append(fits.PrimaryHDU(data=None, header=hdr))
    hdus.append(fits.BinTableHDU(maps, name='MAPS'))
    hdus.writeto(output_fn, overwrite=False)

######################################## Create all-sky map for reproducing the C_ell plots ####################################################

nside = 512
maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/smoothed/desi_dust_gr_{}_allsky.fits'.format(nside)))
maps_sfd = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/v1/smoothed/v1_sfd_ebv_gr_filled_{}.fits'.format(nside)))
assert np.all(maps['HPXPIXEL']==maps_sfd['HPXPIXEL'])
maps['EBV_SFD_INPAINTED'] = maps_sfd['EBV_SFD']
maps['EBV_DIFF'] = maps['EBV_GR']-maps['EBV_SFD_INPAINTED']
maps['RA'], maps['DEC'] = hp.pix2ang(nside, maps['HPXPIXEL'], nest=False, lonlat=True)
maps = maps[['HPXPIXEL', 'RA', 'DEC', 'N_STAR', 'EBV_GR', 'EBV_GR_ERR', 'EBV_SFD_INPAINTED', 'FILL_FRAC']]

# use int32 and float32 to reduce file size
for col in maps.colnames:
    if col in ['HPXPIXEL', 'RA', 'DEC']:
        continue
    if maps[col].dtype=='>i8' or maps[col].dtype=='int64':
        maps[col] = np.array(maps[col], dtype='int32')
    elif maps[col].dtype=='>f8' or maps[col].dtype=='float64':
        maps[col] = np.array(maps[col], dtype='float32')

output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/final_maps/desi_dust_gr_{}_allsky_inpainted.fits'.format(nside)
# maps.write(output_fn, overwrite=False)
# Add header and extension names
hdus = fits.HDUList()
hdr = fits.Header()
hdr['NSIDE'] = nside
hdr['NEST'] = False
hdus.append(fits.PrimaryHDU(data=None, header=hdr))
hdus.append(fits.BinTableHDU(maps, name='MAPS'))
hdus.writeto(output_fn, overwrite=False)


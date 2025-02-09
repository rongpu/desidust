# Created North+South combined catalog and maps
# Correct for the North/South zero point offset

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
import healpy as hp

from sklearn.neighbors import NearestNeighbors
from multiprocessing import Pool

nmad = lambda x: 1.4826*np.median(np.abs(x-np.median(x)))

n_processes = 128

# DECam (airmass 1.3): https://www.legacysurvey.org/dr9/catalogs/#galactic-extinction-coefficients
# BASS+MzLS (airmass 1.1): https://desi.lbl.gov/trac/wiki/ImagingStandardBandpass#SFDExtinctionCoefficients
r = {'g_south': 3.214, 'r_south': 2.165, 'z_south': 1.211, 'g_north': 3.258, 'r_north': 2.176, 'z_north': 1.199}

# From zero_point_north_vs_south-map_based.ipynb
north_gr_offset = -9.5 / 1e3


def hlmean(data, maxpairs=1e8, random_seed=None, verbose=True):
    '''
    Hodges-Lehmann estimator.
    '''

    import itertools

    maxpairs = int(maxpairs)
    ndata = len(data)

    if ndata==0:
        if verbose: print('H-L mean: empty array!!!')
        return None
    if ndata==1:
        return data[0]
    if ndata*(ndata-1)/2 <= maxpairs:
        # only non-identical indices are included
        idx1, idx2 = np.array(list(itertools.combinations(np.arange(ndata), 2))).transpose()
        pairmean1 = (data[idx1]+data[idx2])/2.
        # the pairwise mean of identical indices
        pairmean2 = data
        pairmean = np.concatenate([pairmean1, pairmean2])
        hlmean_value = np.median(pairmean)
    else:
        if verbose: print('Too many pairs; only computing {} pairs'.format(maxpairs))
        if random_seed is not None:
            np.random.seed(random_seed)
        idx1, idx2 = np.random.choice(ndata, size=(maxpairs, 2)).transpose()
        pairmean = (data[idx1]+data[idx2])/2.
        hlmean_value = np.median(pairmean)

    return(hlmean_value)


def get_stats_in_pixel(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    # hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    hp_table['delta_gr_wmean'] = 0.
    hp_table['delta_gr_err'] = 0.
    hp_table['n_star'] = 0
    hp_table['EBV_SFD'] = 0.
    # Additional columns for diagnostics only
    hp_table['delta_gr_mean'] = 0.
    hp_table['delta_gr_median'] = 0.
    hp_table['delta_gr_hlmean'] = 0.

    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        v = cat['dgr'][idx].copy()
        ivar = 1/cat['gr_err'][idx]**2
        hp_table['delta_gr_wmean'][index] = np.sum(v*ivar)/np.sum(ivar)
        hp_table['delta_gr_err'][index] = 1/np.sqrt(np.sum(ivar))
        hp_table['n_star'][index] = len(v)
        hp_table['EBV_SFD'][index] = np.sum(cat['EBV_SFD'][idx]*ivar)/np.sum(ivar)
        hp_table['delta_gr_mean'][index] = np.mean(v)
        hp_table['delta_gr_hlmean'][index] = hlmean(v)
        hp_table['delta_gr_median'][index] = np.median(v)

    return hp_table


catn = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/gr_corrected/gr_north.fits'))
catn1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/gr_corrected/gr_north_predicted_error.fits'))
assert len(catn)==len(catn1) and np.all(catn['TARGETID_DR9']==catn1['TARGETID_DR9'])
catn = hstack([catn, catn1[['gr_err']]])
cats = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/gr_corrected/gr_south.fits'))
cats1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/gr_corrected/gr_south_predicted_error.fits'))
assert len(cats)==len(cats1) and np.all(cats['TARGETID_DR9']==cats1['TARGETID_DR9'])
cats = hstack([cats, cats1[['gr_err']]])

catn['dgr_north'] = catn['dgr'].copy()  # reddening measurement with Northern imaging

# Convert to Southern reddening measurement: remove offset and renormalize
catn['dgr'] = (catn['dgr_north']-north_gr_offset)/(r['g_north']-r['r_north'])*(r['g_south']-r['r_south'])

cats['dgr_north'] = np.nan

mask = catn['DEC']>32.375
catn = catn[mask]

mask = (cats['DEC']<=32.375) | (cats['RA']<70) | (cats['RA']>300)
cats = cats[mask]

cat = vstack([catn, cats], join_type='exact')
print(len(cat))

# # Apply the extra cuts
# mask = cat['extra_cuts'].copy()
# print(np.sum(~mask), np.sum(~mask)/len(mask))
# cat = cat[mask]
# print(len(cat))

mask = cat['PARALLAX']<1.
print('PARALLAX', np.sum(mask)/len(mask))
cat = cat[mask]
print(len(cat))

mask = cat['SN_B']>5
print('SN_B', np.sum(mask)/len(mask))
cat = cat[mask]
print(len(cat))

mask = cat['model_gr_std']<0.03
print('model_gr_std', np.sum(mask)/len(mask))
cat = cat[mask]
print(len(cat))

# Remove duplicates
print(len(cat), len(np.unique(cat['TARGETID_DR9'])))
_, idx = np.unique(cat['TARGETID_DR9'], return_index=True)
cat = cat[idx]
print(len(cat), len(np.unique(cat['TARGETID_DR9'])))

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/gr_corrected/gr_combined.fits', overwrite=False)

############################################################################################################################

# cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/gr_corrected/gr_combined.fits'))

for nside in [64, 128, 256, 512]:

    npix = hp.nside2npix(nside)

    pix_allobj = hp.pixelfunc.ang2pix(nside, cat['TARGET_RA'], cat['TARGET_DEC'], lonlat=True)
    pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

    pixcnts = np.insert(pixcnts, 0, 0)
    pixcnts = np.cumsum(pixcnts)

    pixorder = np.argsort(pix_allobj)

    # split among the Cori processors
    pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_stats_in_pixel, pix_idx_split)

    maps = vstack(res)
    maps.sort('HPXPIXEL')

    maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/dgr_map_combined_{}.fits'.format(nside), overwrite=False)

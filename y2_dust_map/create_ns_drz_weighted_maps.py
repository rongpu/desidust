# Create separate North and South weighted maps (e.g., for determining the North/South relative zero point offset)

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

    hp_table['delta_rz_wmean'] = 0.
    hp_table['delta_rz_err'] = 0.
    hp_table['n_star'] = 0
    hp_table['EBV_SFD'] = 0.
    # Additional columns for diagnostics only
    hp_table['delta_rz_mean'] = 0.
    hp_table['delta_rz_median'] = 0.
    hp_table['delta_rz_hlmean'] = 0.

    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        v = cat['drz'][idx].copy()
        ivar = 1/cat['rz_err'][idx]**2
        hp_table['delta_rz_wmean'][index] = np.sum(v*ivar)/np.sum(ivar)
        hp_table['delta_rz_err'][index] = 1/np.sqrt(np.sum(ivar))
        hp_table['n_star'][index] = len(v)
        hp_table['EBV_SFD'][index] = np.sum(cat['EBV_SFD'][idx]*ivar)/np.sum(ivar)
        hp_table['delta_rz_mean'][index] = np.mean(v)
        hp_table['delta_rz_hlmean'][index] = hlmean(v)
        hp_table['delta_rz_median'][index] = np.median(v)

    return hp_table


for field in ['north', 'south']:

    cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/rz_corrected/rz_{}.fits'.format(field)))
    cat1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/rz_corrected/rz_{}_predicted_error.fits'.format(field)))
    assert len(cat)==len(cat1) and np.all(cat['TARGETID_DR9']==cat1['TARGETID_DR9'])
    cat = hstack([cat, cat1[['rz_err']]])

    mask = cat['PARALLAX']<1.
    print('PARALLAX', np.sum(mask)/len(mask))
    cat = cat[mask]
    print(len(cat))

    mask = cat['SN_B']>5
    print('SN_B', np.sum(mask)/len(mask))
    cat = cat[mask]
    print(len(cat))

    mask = cat['model_rz_std']<0.035
    print('model_rz_std', np.sum(mask)/len(mask))
    cat = cat[mask]
    print(len(cat))

    # Remove duplicates
    print(len(cat), len(np.unique(cat['TARGETID_DR9'])))
    _, idx = np.unique(cat['TARGETID_DR9'], return_index=True)
    cat = cat[idx]
    print(len(cat), len(np.unique(cat['TARGETID_DR9'])))

    ############################################################################################################################

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

        maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_{}_{}.fits'.format(field, nside), overwrite=False)




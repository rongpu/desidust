# North and South separately
# Correct for the per-star_type offsets (relative to bright_mws)
# Create delta_r-z catalog with correction based on the reference catalog (and independent from SFD)
# The extra quality cuts –– the parallax and SN_B –– are included as flags but not applied
# Create unweighted maps

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

# rz_offset_dict = {'dark_std': -1.7, 'bright_std': -2.0, 'bright_mws': 0.0, 'backup_std': -2.7, 'backup_mws': -3.9}  # per_star_type_rz_offsets.ipynb (in mmag)


def get_nn_mean(data):
    tt = Table()
    tt['model_rz_offset'], tt['model_rz_nmad'], tt['model_rz_std'] = np.zeros((3, len(data)))
    for index in range(len(data)):
        idx = nbrs.kneighbors(data[[index]])[1][0]
        v = cat_ref['drz_ref'][idx].copy()  # offset in delta_r-z
        # 5-sigma clipping
        sigma = nmad(v)
        mask = (v<(np.mean(v)-5*sigma)) | (v>(np.mean(v)+5*sigma))
        v = np.clip(v, np.mean(v)-5*sigma, np.mean(v)+5*sigma)
        if np.sum(mask)/len(mask)>0.05:
            print('More than 5% outliers: {}'.format(np.sum(mask)),
                  data[index]*np.array([x1_sigma, x2_sigma, x3_sigma]))
        tt['model_rz_offset'][index] = np.mean(v)
        tt['model_rz_nmad'][index] = sigma
        tt['model_rz_std'][index] = np.std(v)

    return tt


def get_nn_mean_wrapper(data, n_processes=128):
    idx_list = np.array_split(np.arange(len(data)), n_processes)
    data_list = []
    for index in range(len(idx_list)):
        data_list.append(data[idx_list[index]])
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_nn_mean, data_list)
    tt = vstack(res, join_type='exact')

    return tt


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

    hp_table['delta_rz_mean'] = 0.
    hp_table['delta_rz_median'] = 0.
    hp_table['delta_rz_hlmean'] = 0.
    hp_table['n_star'] = 0
    hp_table['EBV_SFD'] = 0.

    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        v = cat['drz'][idx].copy()
        hp_table['delta_rz_mean'][index] = np.mean(v)
        hp_table['delta_rz_hlmean'][index] = hlmean(v)
        hp_table['delta_rz_median'][index] = np.median(v)
        hp_table['n_star'][index] = len(v)
        hp_table['EBV_SFD'][index] = np.mean(cat['EBV_SFD'][idx])

    return hp_table


for field in ['south', 'north']:

    ################################################### Create catalogs ###################################################

    print('\n##################################', field, '##################################')

    cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_combined_{}.fits'.format(field)))
    cat1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_combined_gaia_ls_photom_{}.fits'.format(field)))
    assert len(cat)==len(cat1) and np.all(cat['TARGETID']==cat1['TARGETID'])
    cat1.remove_column('TARGETID')
    cat = hstack([cat, cat1])

    cat.rename_column('EBV', 'EBV_SFD')

    # No extinction correction
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cat['gmag'] = 22.5 - 2.5*np.log10(cat['FLUX_G'])
        cat['rmag'] = 22.5 - 2.5*np.log10(cat['FLUX_R'])
        cat['zmag'] = 22.5 - 2.5*np.log10(cat['FLUX_Z'])

    # Use the SFD-based reddened delta_r-z value for the initial correction and for anchoring the zero point
    model_rz_sfd = np.zeros(len(cat))
    model_rz_sfd[cat['PHOTSYS_PHOTOM']=='N'] = (cat['MODEL_RMAG_N_REDDENED']-cat['MODEL_ZMAG_N_REDDENED'])[cat['PHOTSYS_PHOTOM']=='N']
    model_rz_sfd[cat['PHOTSYS_PHOTOM']=='S'] = (cat['MODEL_RMAG_S_REDDENED']-cat['MODEL_ZMAG_S_REDDENED'])[cat['PHOTSYS_PHOTOM']=='S']
    cat['drz_sfd_raw'] = (cat['rmag']-cat['zmag']) - model_rz_sfd
    assert np.sum(model_rz_sfd==0)==0

    # Use the un-reddened delta_r-z value for the (iterative) correction
    model_rz = np.zeros(len(cat))
    model_rz[cat['PHOTSYS_PHOTOM']=='N'] = (cat['MODEL_RMAG_N']-cat['MODEL_ZMAG_N'])[cat['PHOTSYS_PHOTOM']=='N']
    model_rz[cat['PHOTSYS_PHOTOM']=='S'] = (cat['MODEL_RMAG_S']-cat['MODEL_ZMAG_S'])[cat['PHOTSYS_PHOTOM']=='S']
    cat['drz_raw'] = (cat['rmag']-cat['zmag']) - model_rz
    assert np.sum(model_rz==0)==0

    # extinction corrected colors
    cat['g-r'] = (cat['gmag'] - 3.214*cat['EBV_SFD']) - (cat['rmag']-2.165*cat['EBV_SFD'])
    cat['r-z'] = (cat['rmag'] - 2.165*cat['EBV_SFD']) - (cat['zmag']-1.211*cat['EBV_SFD'])

    cat['EFFTIME_LRG'] = cat['TSNR2_LRG'] * 12.15

    ###################################### Initial quality cuts ######################################

    mask = cat['PHOT_G_MEAN_MAG']!=0
    print('PHOT_G_MEAN_MAG', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['PARALLAX_ERROR']!=0
    print('PARALLAX_ERROR', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    # Remove saturated objects
    mask = (cat['ANYMASK_R']==0) & (cat['ANYMASK_Z']==0)
    print('ANYMASK_R/Z', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['ZWARN']==0
    print('ZWARN', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['DELTACHI2']>100
    print('DELTACHI2', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['SPECTYPE']=='STAR'
    print('SPECTYPE', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['SUBTYPE']!='WD'
    print('SUBTYPE', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['RVS_WARN']==0
    print('RVS_WARN', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    # Remove galaxies
    mask = cat['Z']<0.002  # There are no galaxies at z<0.002 based on VI
    print('Z', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    mask = cat['EFFTIME_LRG']>30
    print('EFFTIME_LRG', np.sum(~mask), len(mask), np.sum(~mask)/len(mask))
    cat = cat[mask]

    # Keep unique objects –– choose the highest BLUE_SN object
    print(len(cat), len(np.unique(cat['TARGETID_DR9'])))
    cat.sort('SN_B', reverse=True)
    _, idx_keep = np.unique(cat['TARGETID_DR9'], return_index=True)
    cat = cat[idx_keep]
    print(len(cat), len(np.unique(cat['TARGETID_DR9'])))

    if field=='north':
        mask = cat['DEC']>30
        print('North Dec cut', np.sum(~mask), len(mask), np.sum(~mask)/len(mask))
        cat = cat[mask]

    ################################################################
    mask_extra = cat['PARALLAX']<1.
    print('PARALLAX', np.sum(~mask_extra), np.sum(~mask_extra)/len(mask_extra))
    mask_extra &= cat['SN_B']>8
    print('+ SN_B', np.sum(~mask_extra), np.sum(~mask_extra)/len(mask_extra))
    cat['extra_cuts'] = mask_extra.copy()
    ################################################################

    print(len(cat))

    print('More quality cuts')

    mask_quality = (cat['NOBS_G']>0) & (cat['NOBS_R']>0) & (cat['NOBS_Z']>0)
    print(np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    mask_quality &= (cat['FLUX_IVAR_G']>0) & (cat['FLUX_IVAR_R']>0) & (cat['FLUX_IVAR_Z']>0)
    print(np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    mask_quality &= (cat['FRACFLUX_R']<0.01) & (cat['FRACFLUX_Z']<0.01)
    print(np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    mask_quality &= (cat['FRACMASKED_R']<0.6) & (cat['FRACMASKED_Z']<0.6)
    print(np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    mask_quality &= (cat['FIBERTOTFLUX_R']>0) & (cat['FIBERTOTFLUX_Z']>0)
    print('FIBERTOTFLUX>0', np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    mask_quality &= (cat['FIBERFLUX_R']/cat['FIBERTOTFLUX_R']>0.99) & (cat['FIBERFLUX_Z']/cat['FIBERTOTFLUX_Z']>0.99)
    print(np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    mask_quality &= (cat['PHOT_G_MEAN_MAG']!=0) & (cat['PHOT_BP_MEAN_MAG']!=0) & (cat['PHOT_RP_MEAN_MAG']!=0)
    print(np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
    cat = cat[mask_quality]
    print(len(cat))

    mask = np.abs(cat['ls_mag_r']-cat['rmag'])<0.1
    mask &= np.abs(cat['ls_mag_z']-cat['zmag'])<0.1
    print('Gaia outliers', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]
    print(len(cat))

    print('Stellar parameter cuts')
    mask_good = (cat['LOGG']>3.5) & (cat['LOGG']<5.2)
    print(np.sum(~mask_good), np.sum(~mask_good)/len(mask_good))
    mask_good &= (cat['TEFF']>5000) & (cat['TEFF']<6500)
    print(np.sum(~mask_good), np.sum(~mask_good)/len(mask_good))
    mask_good &= (cat['FEH']>-4.0) & (cat['FEH']<0.5)
    print(np.sum(~mask_good), np.sum(~mask_good)/len(mask_good))
    cat = cat[mask_good]
    print(len(cat))

    cat_ref = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/rz_corrected/rz_reference_{}.fits'.format(field)))

    # Apply random downselection on the reference catalog but not the catalog used to generate the EBV maps
    mask = cat_ref['downselect'].copy()
    cat_ref = cat_ref[mask]
    cat_ref.remove_column('downselect')
    print(len(cat_ref))

    x1, x2, x3 = cat_ref['LOGG'].copy(), cat_ref['FEH'].copy(), cat_ref['TEFF'].copy()
    x1_sigma, x2_sigma, x3_sigma = nmad(x1), nmad(x2), nmad(x3)
    x1 /= x1_sigma
    x2 /= x2_sigma
    x3 /= x3_sigma
    X = np.vstack([x1, x2, x3]).T
    nbrs = NearestNeighbors(n_neighbors=400, algorithm='ball_tree').fit(X)

    X_all = np.vstack([cat['LOGG']/x1_sigma, cat['FEH']/x2_sigma, cat['TEFF']/x3_sigma]).T
    tt = get_nn_mean_wrapper(X_all, n_processes=n_processes)
    if 'model_rz_offset' in cat.colnames:
        cat.remove_columns(['model_rz_offset', 'model_rz_nmad', 'model_rz_std'])
    cat = hstack([cat, tt], join_type='exact')

    # offset-corrected delta_r-z
    cat['drz'] = cat['drz_raw'] - cat['model_rz_offset']  # To get the corrected model_rz: model_rz_corr = model_rz + model_rz_offset

    # # correct for per-star_type offsets
    # drz_new = np.full(len(cat), -99.)
    # for star_type in ['dark_std', 'bright_std', 'bright_mws', 'backup_std', 'backup_mws']:
    #     # beware of objects that are in multiple star types
    #     mask = cat[star_type].copy()
    #     drz_new[mask] = cat['drz'][mask] - 1e-3*rz_offset_dict[star_type]
    # assert np.all(drz_new!=-99.)
    # cat['drz'] = drz_new

    cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/rz_corrected/rz_{}.fits'.format(field), overwrite=True)

    ################################################### Create maps ###################################################

    # Apply the extra cuts
    mask = cat['extra_cuts'].copy()
    print('Extra cuts', np.sum(~mask), np.sum(~mask)/len(mask))
    cat = cat[mask]
    print(len(cat))

    # Remove duplicates
    print(len(cat), len(np.unique(cat['TARGETID_DR9'])))
    _, idx = np.unique(cat['TARGETID_DR9'], return_index=True)
    cat = cat[idx]
    print(len(cat), len(np.unique(cat['TARGETID_DR9'])))

    for nside in [64, 128, 256]:

        pix_allobj = hp.pixelfunc.ang2pix(nside, cat['TARGET_RA'], cat['TARGET_DEC'], lonlat=True)
        pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)
        pixcnts = np.insert(pixcnts, 0, 0)
        pixcnts = np.cumsum(pixcnts)
        pixorder = np.argsort(pix_allobj)
        pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)
        with Pool(processes=n_processes) as pool:
            res = pool.map(get_stats_in_pixel, pix_idx_split)

        maps = vstack(res)
        maps.sort('HPXPIXEL')

        maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_{}_{}.fits'.format(field, nside), overwrite=True)

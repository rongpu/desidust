# Add Gaia DR3 and DR9 sweep columns
# Assign star type
# Combine Y1 and Y1-2 daily data
# Create North-only and South-only catalogs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from desitarget import targetmask

# Exclude Gaia columns
columns = ['TARGETID', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'FA_TARGET', 'FA_TYPE', 'OBJTYPE', 'SUBPRIORITY', 'OBSCONDITIONS', 'RELEASE', 'BRICKNAME', 'BRICKID', 'BRICK_OBJID', 'MORPHTYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'MASKBITS', 'SERSIC', 'SHAPE_R', 'SHAPE_E1', 'SHAPE_E2', 'PHOTSYS', 'PRIORITY_INIT', 'NUMOBS_INIT', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'PLATE_RA', 'PLATE_DEC', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y', 'RMS_DELTA_Y', 'MEAN_FIBER_RA', 'STD_FIBER_RA', 'MEAN_FIBER_DEC', 'STD_FIBER_DEC', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'VRAD', 'VRAD_ERR', 'VRAD_SKEW', 'VRAD_KURT', 'LOGG', 'TEFF', 'ALPHAFE', 'FEH', 'LOGG_ERR', 'TEFF_ERR', 'ALPHAFE_ERR', 'FEH_ERR', 'VSINI', 'CHISQ_TOT', 'CHISQ_C_TOT', 'CHISQ_B', 'CHISQ_C_B', 'CHISQ_R', 'CHISQ_C_R', 'CHISQ_Z', 'CHISQ_C_Z', 'RVS_WARN', 'SN_B', 'SN_R', 'SN_Z', 'SUCCESS', 'RR_Z', 'RR_SPECTYPE', 'MODEL_GMAG_N_REDDENED', 'MODEL_GMAG_N', 'MODEL_GMAG_S_REDDENED', 'MODEL_GMAG_S', 'MODEL_RMAG_N_REDDENED', 'MODEL_RMAG_N', 'MODEL_RMAG_S_REDDENED', 'MODEL_RMAG_S', 'MODEL_IMAG_S_REDDENED', 'MODEL_IMAG_S', 'MODEL_ZMAG_N_REDDENED', 'MODEL_ZMAG_N', 'MODEL_ZMAG_S_REDDENED', 'MODEL_ZMAG_S', 'program', 'SCND_TARGET', 'Z', 'ZERR', 'ZWARN', 'CHI2', 'NPIXELS', 'SPECTYPE', 'SUBTYPE', 'DELTACHI2', 'TSNR2_GPBDARK', 'TSNR2_ELG', 'TSNR2_GPBBRIGHT', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_GPBBACKUP', 'TSNR2_QSO', 'TSNR2_LRG']

input_paths = [
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_y1.fits',
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_daily.fits',
]

cat_dict = {}

for input_path in input_paths:

    if 'stars_daily.fits' in input_path:
        cat = Table(fitsio.read(input_path, columns=columns+['TILEID']))
    else:
        cat = Table(fitsio.read(input_path, columns=columns))

    # Use Gaia DR3
    gaia_path = input_path.replace('.fits', '-gaia_dr3.fits')
    gaia = Table(fitsio.read(gaia_path))
    mask = gaia['REF_EPOCH']!=0
    print(np.sum(mask), np.sum(mask)/len(mask))
    gaia = gaia[mask]
    cat = cat[mask]
    cat = hstack([cat, gaia])

    # Remove objects from low EFFTIME tiles
    if 'stars_daily.fits' in input_path:
        tiles = Table.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/misc/tiles-daily.csv')
        tiles1 = Table.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/misc/tiles-main.ecsv')
        tiles1 = tiles1[tiles1['IN_DESI']]
        mask = np.in1d(tiles['TILEID'], tiles1['TILEID'])
        tiles = tiles[mask]
        cat = join(cat, tiles[['TILEID', 'EFFTIME_SPEC']], join_type='left')
        mask = cat['EFFTIME_SPEC']>40
        cat = cat[mask]
        cat.remove_columns(['TILEID', 'EFFTIME_SPEC'])

    print(len(cat), len(np.unique(cat['SOURCE_ID'])), len(np.unique(cat['SOURCE_ID']))/len(cat))

    cat_all = cat.copy()

    for field in ['north', 'south']:

        cat = cat_all.copy()

        sweep_path = input_path.replace('.fits', '-dr9_sweep_{}.fits'.format(field))

        sweep = Table(fitsio.read(sweep_path))
        # remove Gaia columns from sweep
        sweep_remove_columns = ['RA', 'DEC', 'RA_IVAR', 'DEC_IVAR', 'REF_CAT', 'REF_ID', 'REF_EPOCH', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_G_MEAN_FLUX_OVER_ERROR', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_FLUX_OVER_ERROR', 'GAIA_PHOT_RP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_FLUX_OVER_ERROR', 'GAIA_ASTROMETRIC_EXCESS_NOISE', 'GAIA_DUPLICATED_SOURCE', 'GAIA_PHOT_BP_RP_EXCESS_FACTOR', 'GAIA_ASTROMETRIC_SIGMA5D_MAX', 'GAIA_ASTROMETRIC_PARAMS_SOLVED', 'PARALLAX', 'PARALLAX_IVAR', 'PMRA', 'PMRA_IVAR', 'PMDEC', 'PMDEC_IVAR']
        sweep.remove_columns(sweep_remove_columns)

        print(len(sweep)==len(np.unique(sweep['TARGETID_DR9'])), len(sweep), len(np.unique(sweep['TARGETID_DR9'])))
        print(len(sweep)==len(np.unique(sweep['TARGETID'])), len(sweep), len(np.unique(sweep['TARGETID'])))

        # Require valid photometry
        mask = (sweep['FLUX_G']>0) & (sweep['FLUX_IVAR_G']>0) & (sweep['FLUX_R']>0) & (sweep['FLUX_IVAR_R']>0)
        print(np.sum(~mask), np.sum(~mask)/len(mask))
        sweep = sweep[mask]

        mask = np.in1d(cat['TARGETID'], sweep['TARGETID'])
        cat = cat[mask]

        # Remove sweep columns from the original catalog
        cols_to_remove = list(np.intersect1d(cat.colnames, sweep.colnames))
        cols_to_remove.remove('TARGETID')
        cat.remove_columns(cols_to_remove)
        print(len(cat))
        cat = join(cat, sweep, keys='TARGETID', join_type='left')
        print(len(cat))

        # ############# Keep unique objects –– choose the highest BLUE_SN object #############
        # cat.sort('SN_B', reverse=True)
        # _, idx_keep = np.unique(cat['TARGETID_DR9'], return_index=True)
        # cat = cat[idx_keep]
        # print(len(cat), len(np.unique(cat['TARGETID_DR9'])))
        # print(len(cat)==len(np.unique(cat['TARGETID_DR9'])), len(cat), len(np.unique(cat['TARGETID_DR9'])))
        # print(len(cat)==len(np.unique(cat['TARGETID'])), len(cat), len(np.unique(cat['TARGETID'])))

        # PHOTSYS_PHOTOM denotes the origin of the photometry; PHOTSYS (from spectroscopic FIBERMAP) is retained but renamed to PHOTSYS_FIBERMAP
        if field=='north':
            cat['PHOTSYS_PHOTOM'] = 'N'
        else:
            cat['PHOTSYS_PHOTOM'] = 'S'
        cat.rename_column('PHOTSYS', 'PHOTSYS_FIBERMAP')

        if 'stars_y1.fits' in input_path:
            cat['y1'] = True
            cat_dict['y1_'+field] = cat.copy()
        else:
            cat['y1'] = False
            cat_dict['daily_'+field] = cat.copy()


for field in ['north', 'south']:

    cat_stack = []
    for key in ['y1', 'daily']:
        cat_stack.append(cat_dict[key+'_'+field].copy())
    cat = vstack(cat_stack, join_type='exact')

    # main dark selection
    cat['dark_std'] = cat['program']=='main dark'
    desi_target_list = ['STD_FAINT', 'STD_BRIGHT']
    mask = np.full(len(cat), False)
    for desi_target in desi_target_list:
        mask |= cat['DESI_TARGET'] & targetmask.desi_mask[desi_target] > 0
    cat['dark_std'] &= mask
    print(np.sum(cat['dark_std']))

    # main bright standards selection
    cat['bright_std'] = cat['program']=='main bright'
    desi_target_list = ['STD_FAINT', 'STD_BRIGHT']
    mask = np.full(len(cat), False)
    for desi_target in desi_target_list:
        mask |= cat['DESI_TARGET'] & targetmask.desi_mask[desi_target] > 0
    cat['bright_std'] &= mask
    print(np.sum(cat['bright_std']))

    # main bright MWS selection
    cat['bright_mws'] = cat['program']=='main bright'
    mws_target_list = ['MWS_MAIN_BLUE', 'MWS_MAIN_RED', 'MWS_BROAD']
    mask = np.full(len(cat), False)
    for mws_target in mws_target_list:
        mask |= cat['MWS_TARGET'] & targetmask.mws_mask[mws_target] > 0
    cat['bright_mws'] &= mask
    print(np.sum(cat['bright_mws']))

    # main backup standards selection
    cat['backup_std'] = cat['program']=='main backup'
    mws_target_list = ['GAIA_STD_FAINT', 'GAIA_STD_BRIGHT']
    mask = np.full(len(cat), False)
    for mws_target in mws_target_list:
        mask |= cat['MWS_TARGET'] & targetmask.mws_mask[mws_target] > 0
    cat['backup_std'] &= mask
    print(np.sum(cat['backup_std']))

    # main backup MWS selection
    cat['backup_mws'] = cat['program']=='main backup'
    mws_target_list = ['BACKUP_BRIGHT', 'BACKUP_FAINT']
    mask = np.full(len(cat), False)
    for mws_target in mws_target_list:
        mask |= cat['MWS_TARGET'] & targetmask.mws_mask[mws_target] > 0
    cat['backup_mws'] &= mask
    print(np.sum(cat['backup_mws']))

    cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_combined_{}.fits'.format(field))

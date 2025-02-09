# source /global/common/software/desi/desi_environment.sh 23.1

# Y1-2 daily reduction
# Combine the coadded spectra, redrock tables and rvspecfit model spectra into a single FITS file

# rongpu@login34:spectra> du -sh *
# 1.2T    daily
# 999G    y1

# The data is archived in HPSS: /home/r/rongpu/desi_dust/spectra.tar


from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool

from desitarget import targetmask

programs = ['main dark', 'main bright', 'main backup']
output_dir = '/pscratch/sd/r/rongpu/ebv/spectra/daily/'


def select_stars(fibermap, program):

    if program=='main dark':
        desi_target_list = ['STD_FAINT', 'STD_BRIGHT']
        mask = np.full(len(fibermap), False)
        for desi_target in desi_target_list:
            mask |= fibermap['DESI_TARGET'] & targetmask.desi_mask[desi_target] > 0

    elif program=='main bright':
        mask = np.full(len(fibermap), False)

        desi_target_list = ['STD_FAINT', 'STD_BRIGHT']
        for desi_target in desi_target_list:
            mask |= fibermap['DESI_TARGET'] & targetmask.desi_mask[desi_target] > 0

        mws_target_list = ['MWS_MAIN_BLUE', 'MWS_MAIN_RED', 'MWS_BROAD']
        for mws_target in mws_target_list:
            mask |= fibermap['MWS_TARGET'] & targetmask.mws_mask[mws_target] > 0

    elif program=='main backup':
        mws_target_list = ['GAIA_STD_FAINT', 'GAIA_STD_BRIGHT', 'BACKUP_BRIGHT', 'BACKUP_FAINT']
        mask = np.full(len(fibermap), False)
        for mws_target in mws_target_list:
            mask |= fibermap['MWS_TARGET'] & targetmask.mws_mask[mws_target] > 0

    else:
        raise ValueError

    return mask


def combine_data(tileid):

    tile_index = np.where(tiles['TILEID']==tileid)[0][0]
    program = tiles['FAFLAVOR'][tile_index].replace('main', 'main ')

    lastnight = tiles['LASTNIGHT'][tile_index]
    tilenight = str(tileid)+'/'+str(lastnight)

    fns = glob.glob('/pscratch/sd/k/koposov/specmod/daily_211013/tiles/cumulative/'+tilenight+'/*.fits')
    fns.sort()

    cat_stack = []

    for fn in fns:

        # remove duplicate tmp files, e.g., /dvs_ro/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/4400/20230508/coadd-0-4400-thru20230421_tmp55110.fits
        if '_tmp' in fn:
            continue

        fn1 = os.path.join('/dvs_ro/cfs/cdirs/desicollab/science/mws/redux/daily/rv_output/231013/tiles/cumulative', tilenight, os.path.basename(fn).replace('specmod_coadd', 'rvtab_coadd'))
        fn2 = os.path.join('/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative', tilenight, os.path.basename(fn).replace('specmod_coadd', 'coadd'))  # daily coadd file
        fn_rr = os.path.join('/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative', tilenight, os.path.basename(fn).replace('specmod_coadd', 'redrock'))

        fibermap = Table(fitsio.read(fn, ext='FIBERMAP'))
        mask = select_stars(fibermap, program)
        fibermap = fibermap[mask]
        if len(fibermap)==0:
            return None

        fn_output = os.path.join(output_dir, tilenight, os.path.basename(fn).replace('specmod_coadd', 'rvspecfit'))
        if not os.path.exists(os.path.dirname(fn_output)):
            try:  # in case another process is also creating the directory
                os.makedirs(os.path.dirname(fn_output))
            except:
                pass
        hdul = fitsio.FITS(fn_output, mode='rw', clobber=True)
        hdul.write(None)

        fn_input = fn2
        tmp = Table(fitsio.read(fn_input, ext='FIBERMAP'))
        idx = np.where(np.in1d(tmp['TARGETID'], fibermap['TARGETID']))[0]
        assert len(fibermap)==len(idx) and np.all(fibermap['TARGETID']==tmp['TARGETID'][idx])
        for ext in ['FIBERMAP', 'B_WAVELENGTH', 'B_FLUX', 'B_IVAR', 'B_MASK', 'R_WAVELENGTH', 'R_FLUX', 'R_IVAR', 'R_MASK', 'Z_WAVELENGTH', 'Z_FLUX', 'Z_IVAR', 'Z_MASK', 'SCORES']:
            hdr = fitsio.read_header(fn_input, ext=ext)
            data = fitsio.read(fn_input, ext=ext)
            if '_WAVELENGTH' not in ext:
                data = data[idx]
            hdul.write(data, header=hdr, extname=ext)

        fn_input = fn_rr
        tmp = Table(fitsio.read(fn_input, ext='FIBERMAP'))
        idx = np.where(np.in1d(tmp['TARGETID'], fibermap['TARGETID']))[0]
        assert len(fibermap)==len(idx) and np.all(fibermap['TARGETID']==tmp['TARGETID'][idx])
        for ext in ['REDSHIFTS']:
            hdr = fitsio.read_header(fn_input, ext=ext)
            data = fitsio.read(fn_input, ext=ext)
            data = data[idx]
            hdul.write(data, header=hdr, extname=ext)

        fn_input = fn1
        tmp = Table(fitsio.read(fn_input, ext='FIBERMAP'))
        idx = np.where(np.in1d(tmp['TARGETID'], fibermap['TARGETID']))[0]
        assert len(fibermap)==len(idx) and np.all(fibermap['TARGETID']==tmp['TARGETID'][idx])
        for ext in ['RVTAB']:
            hdr = fitsio.read_header(fn_input, ext=ext)
            data = fitsio.read(fn_input, ext=ext)
            data = data[idx]
            hdul.write(data, header=hdr, extname=ext)

        fn_input = fn
        tmp = Table(fitsio.read(fn_input, ext='FIBERMAP'))
        idx = np.where(np.in1d(tmp['TARGETID'], fibermap['TARGETID']))[0]
        assert len(fibermap)==len(idx) and np.all(fibermap['TARGETID']==tmp['TARGETID'][idx])
        for ext in ['WAVELENGTH', 'SPECMOD']:
            hdr = fitsio.read_header(fn_input, ext=ext)
            data = fitsio.read(fn_input, ext=ext)
            if ext!='WAVELENGTH':
                data = data[idx]
            hdul.write(data, header=hdr, extname=ext)

        hdul.close()

    return None


# tiles = Table.read('/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/tiles-daily.csv')
# tiles1 = Table.read('/dvs_ro/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv')
tiles = Table.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/misc/tiles-daily.csv')
tiles1 = Table.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/misc/tiles-main.ecsv')
tiles1 = tiles1[tiles1['IN_DESI']]  # Exclude tertiary tiles, calibration tiles, etc.
mask = np.in1d(tiles['TILEID'], tiles1['TILEID'])
tiles = tiles[mask]
print(len(tiles))
mask = tiles['EFFTIME_SPEC']>0
tiles = tiles[mask]
print(len(tiles))

tileid_list = []
for index in range(len(tiles)):
    tileid, lastnight = tiles['TILEID'][index], tiles['LASTNIGHT'][index]
    fns = glob.glob('/pscratch/sd/k/koposov/specmod/daily_211013/tiles/cumulative/{}/{}/*.fits'.format(tileid, lastnight))
    if len(fns)>0:
        tileid_list.append(tileid)

n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(combine_data, tileid_list, chunksize=1)


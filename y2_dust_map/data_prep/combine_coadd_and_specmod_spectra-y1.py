# source /global/common/software/desi/desi_environment.sh 23.1

# Y1 iron reduction
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
output_dir = '/pscratch/sd/r/rongpu/ebv/spectra/y1/'


def select_stars(fibermap):

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

    return mask


def combine_data(fn):

    a, b, hpix = fn.split('/')[-5], fn.split('/')[-4], int(fn.split('/')[-2])
    fn1 = '/dvs_ro/cfs/cdirs/desicollab/science/mws/redux/iron/rv_output/230825/healpix/{}/{}/{}/{}/rvtab_coadd-{}-{}-{}.fits'.format(a, b, hpix//100, hpix, a, b, hpix)
    fn2 = fn.replace('/pscratch/sd/k/koposov/specmod/iron_230825/', '/dvs_ro/cfs/cdirs/desi/spectro/redux/iron/healpix/').replace('specmod_', '')  # iron coadd file
    fn_rr = fn.replace('/pscratch/sd/k/koposov/specmod/iron_230825/', '/dvs_ro/cfs/cdirs/desi/spectro/redux/iron/healpix/').replace('specmod_coadd', 'redrock')

    fibermap = Table(fitsio.read(fn, ext='FIBERMAP'))
    mask = select_stars(fibermap)
    fibermap = fibermap[mask]
    if len(fibermap)==0:
        return None

    fn_output = os.path.join(output_dir, '{}/{}/{}/{}/rvspecfit-{}-{}-{}.fits'.format(a, b, hpix//100, hpix, a, b, hpix))
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


for program in programs:

    a, b = program.split()[0], program.split()[1]
    fns = glob.glob('/pscratch/sd/k/koposov/specmod/iron_230825/{}/{}/*/*/*.fits'.format(a, b))
    print(program, len(fns))
    n_process = 256
    with Pool(processes=n_process) as pool:
        res = pool.map(combine_data, fns, chunksize=1)

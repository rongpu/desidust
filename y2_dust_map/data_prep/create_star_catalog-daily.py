# Y1-2 daily reduction
# Compute no-reddening and (SFD-)reddened magnitudes in BASS/MzLS and DECam filters
# iron tag: "source /global/common/software/desi/desi_environment.sh 23.1"
# desispec/0.56.5
# speclite/0.16

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool

import speclite.filters
from desispec.magnitude import compute_ab_mag
from desiutil.dust import dust_transmission

programs = ['main dark', 'main bright', 'main backup']


def load_legacy_survey_filter(band, photsys):
    """
    Uses speclite.filters to load the filter transmission
    Returns speclite.filters.FilterResponse object

    Args:
        band: filter pass-band in "G","R","Z","W1","W2"
        photsys: "N" or "S" for North (BASS+MzLS) or South (CTIO/DECam)
    """
    filternamemap=None
    if band[0].upper()=="W":  # it's WISE
        filternamemap = "wise2010-{}".format(band.upper())
    elif band.upper() in ['G', 'R', 'I', 'Z']:
        if photsys=="N":
            if band.upper() in ["G", "R"]:
                filternamemap="BASS-{}".format(band.lower())
            else:
                filternamemap="MzLS-z"
        elif photsys=="S":
            ##################################################
            filternamemap="decamDR1-{}".format(band.lower())
            ##################################################
        else:
            raise ValueError("unknown photsys '{}', known ones are 'N' and 'S'".format(photsys))
    else:
        raise ValueError("unknown band '{}', known ones are 'G','R', 'I', 'Z','W1' and 'W2'".format(photsys))

    filter_response = speclite.filters.load_filter(filternamemap)
    return filter_response


filter_curves = dict()
for band in ['G', 'R', 'I', 'Z']:
    for photsys in ['N', 'S']:
        filtername = band + photsys
        if filtername=='IN':  # No i band in North
            continue
        filter_curves[filtername] = load_legacy_survey_filter(band=band, photsys=photsys)
# for band in ['W1', 'W2']:
#     filter_curves[band] = load_legacy_survey_filter(band=band, photsys=None)


def get_data(fn):

    cat1 = Table(fitsio.read(fn, ext='FIBERMAP'))
    cat2 = Table(fitsio.read(fn, ext='RVTAB'))
    cat3 = Table(fitsio.read(fn, ext='SCORES'))
    cat4 = Table(fitsio.read(fn, ext='REDSHIFTS'))
    cat2.remove_columns(np.intersect1d(cat1.colnames, cat2.colnames))
    cat3.remove_column('TARGETID')
    cat4.remove_column('TARGETID')
    cat = hstack([cat1, cat2, cat3, cat4], join_type='exact')
    cat['HPXPIXEL'] = int(fn.split('/')[-2])

    stdwave = fitsio.read(fn, ext='WAVELENGTH')
    # print(stdwave.shape)
    models = fitsio.read(fn, ext='SPECMOD')
    # print(models.shape)

    for band in ['G', 'R', 'I', 'Z']:
        for photsys in ['N', 'S']:
            if band+photsys=='IN':  # No i band in North
                continue
            cat['MODEL_{}MAG_{}_REDDENED'.format(band, photsys)] = 0.
            cat['MODEL_{}MAG_{}'.format(band, photsys)] = 0.

    for index in range(len(cat)):

        ebv = cat['EBV'][index]
        # photsys = cat['PHOTSYS'][index]

        model_no_reddening = models[index].copy()
        model = model_no_reddening * dust_transmission(stdwave, ebv)  # apply dereddening

        for band in ['G', 'R', 'I', 'Z']:
            for photsys in ['N', 'S']:
                filtername = band + photsys
                if filtername=='IN':  # No i band in North
                    continue
                filt_ww, filt_tt = filter_curves[filtername].wavelength.copy(), filter_curves[filtername].response.copy()
                # Reddened model magnitude using Julien's script
                cat['MODEL_{}MAG_{}_REDDENED'.format(band, photsys)][index] = compute_ab_mag(stdwave, model, filt_ww, filt_tt)
                # Zero-reddening model magnitude using Julien's script
                cat['MODEL_{}MAG_{}'.format(band, photsys)][index] = compute_ab_mag(stdwave, model_no_reddening, filt_ww, filt_tt)

    # Example /pscratch/sd/r/rongpu/ebv/spectra_new/daily/2770/20230617/rvspecfit-6-2770-thru20230617.fits
    tileid = int(fn.split('/')[-3])
    assert fn.split('/')[-4]=='daily'
    tile_index = np.where(tiles['TILEID']==tileid)[0][0]
    program = tiles['FAFLAVOR'][tile_index].replace('main', 'main ')
    cat['program'] = program

    return cat


tiles = Table.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/misc/tiles-daily.csv')

fns = glob.glob('/pscratch/sd/r/rongpu/ebv/spectra/daily/*/*/*.fits')
fns.sort()
print(len(fns))

n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(get_data, fns, chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)
cat = vstack(res).filled(0)

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_daily.fits')


from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import sys, os, time, argparse, glob
import fitsio
import gc
from astropy.table import Table, vstack, join
from multiprocessing import Pool
import healpy as hp

from desitarget.targets import decode_targetid, encode_targetid

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, scatter_plot

time_start = time.perf_counter()


def do_stuff(hpix):

    mask = cat['hpix']==hpix
    cat2 = cat[['RA', 'DEC', 'id']][mask].copy()

    gaia_fn = '/dvs_ro/cfs/cdirs/cosmo/data/gaia/dr3/healpix/healpix-{}.fits'.format(str(hpix).zfill(5))
    cat1 = Table(fitsio.read(gaia_fn))

    ra1 = np.array(cat1['RA'])
    dec1 = np.array(cat1['DEC'])

    ra2 = np.array(cat2['RA'])
    dec2 = np.array(cat2['DEC'])

    search_radius = 0.1
    idx1, idx2, d2d, d_ra, d_dec = match_coord(ra1, dec1, ra2, dec2, search_radius=search_radius, plot_q=False, keep_all_pairs=True, verbose=False)

    cat1 = cat1[idx1].copy()
    cat1['id'] = cat2['id'][idx2].copy()
    cat1['d2d'] = d2d

    return cat1


input_paths = [
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_y1.fits',
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_daily.fits',
]

# gaia_fns = sorted(glob.glob('/global/cfs/cdirs/cosmo/data/gaia/dr3/healpix/healpix-*.fits'))

for input_path in input_paths:

    output_path = input_path.replace('.fits', '-gaia_dr3.fits').replace('/dvs_ro/', '/global/')

    cat = Table(fitsio.read(input_path))
    print(len(cat), len(np.unique(cat['TARGETID'])))

    # The reference epoch for Gaia DR3 (both Gaia EDR3 and the full Gaia DR3) is 2016.0
    # https://www.cosmos.esa.int/web/gaia/dr3
    cat['RA'] = cat['TARGET_RA'] - (cat['REF_EPOCH'] - 2016.0) * cat['PMRA'] * 1e-3/3600 / np.cos(np.radians(cat['TARGET_DEC']))
    cat['RA'] = (cat['RA'] + 360)%360  # Wrap around
    cat['DEC'] = cat['TARGET_DEC'] - (cat['REF_EPOCH'] - 2016.0) * cat['PMDEC'] * 1e-3/3600

    cat['id'] = np.arange(len(cat))
    cat['hpix'] = hp.ang2pix(32, cat['RA'], cat['DEC'], nest=True, lonlat=True)
    hpix_list = np.unique(cat['hpix'])

    n_process = 256
    with Pool(processes=n_process) as pool:
        res = pool.map(do_stuff, hpix_list)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    gaia = vstack(res)

    gaia = join(cat[['id']], gaia, keys='id', join_type='left').filled(0)
    gaia.sort('id')
    assert np.all(gaia['id']==cat['id'])
    gaia.remove_column('id')

    gaia.write(output_path)


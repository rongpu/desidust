from __future__ import division, print_function
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys, os, time, argparse, glob
import fitsio
import gc
from astropy.table import Table, vstack
from multiprocessing import Pool

from desitarget.targets import decode_targetid, encode_targetid

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, scatter_plot

time_start = time.perf_counter()


def do_stuff(cat1_index):

    # Load sweep catalog

    cat1_path = cat1_paths[cat1_index]
    filename = cat1_path[-26:-5]

    # Match only overlapping regions to reduce computation time
    # Area of the brick from the filename
    # # Reduced sweep files:
    # brick = cat1_path[-28:-13]
    # Original sweep files:
    brick = cat1_path[-20:-5]
    ra1min = float(brick[0:3])
    ra1max = float(brick[8:11])
    dec1min = float(brick[4:7])
    if brick[3]=='m':
        dec1min = -dec1min
    dec1max = float(brick[-3:])
    if brick[-4]=='m':
        dec1max = -dec1max
    mask = (ra2full<ra1max+1/60.) & (ra2full>ra1min-1/60.) & (dec2full<dec1max+1/360.) & (dec2full>dec1min-1/360.)
    if np.sum(mask)==0:
        # print('0 matches')
        return None

    # cat2_idx keeps track of cat2 original index
    cat2_idx = np.arange(len(cat2))
    # keep only cat2 objects in the overlapping region
    ra2 = ra2full[mask]
    dec2 = dec2full[mask]
    cat2_idx = cat2_idx[mask]
    # print('%d out of %d objects in cat2 are in the overlapping region'%(np.sum(mask), len(mask)))

    cat1 = Table(fitsio.read(cat1_path, ext=1, columns=['TYPE', 'GAIA_PHOT_G_MEAN_MAG']))
    # # Remove "DUP" objects
    # mask = (cat1['TYPE']!='DUP') & (cat1['TYPE']!='DUP ')
    # PSF objects only
    mask = cat1['TYPE']=='PSF'
    # Remove objects not in Gaia
    mask &= cat1['GAIA_PHOT_G_MEAN_MAG']!=0
    idx = np.where(mask)[0]
    if len(idx)==0:
        return None
    cat1 = Table(fitsio.read(cat1_path, ext=1, rows=idx))

    ra1 = np.array(cat1['RA'])
    dec1 = np.array(cat1['DEC'])

    idx1, idx2, d2d, d_ra, d_dec = match_coord(ra1, dec1, ra2, dec2, search_radius=search_radius, plot_q=False, keep_all_pairs=True, verbose=False)
    if np.sum(d2d>0.01)>0:
        print('%d - '%cat1_index + filename)
        print('!!! {} objects with >0.01 arcsec match !!'.format(np.sum(d2d>0.01)), cat1_path)
        print('max(d_ra), max(d_dec)', np.max(np.abs(d_ra)), np.max(np.abs(d_dec)))

        # markersize = np.max([0.01, np.min([5, 0.2*100000/len(d_ra)])])
        # axis = [-search_radius*1.05, search_radius*1.05,
        #         -search_radius*1.05, search_radius*1.05]
        # ax = scatter_plot(d_ra, d_dec, markersize=markersize, alpha=0.4, axis=axis, show=False)
        # # plt.savefig(os.path.join(plot_path, '{}_{}.png'.format(cat2_index, brick)))
        # plt.show()

    if len(idx1)==0:
        return None

    idx2_original = cat2_idx[idx2]

    cat1 = cat1[idx1].copy()
    cat1['idx_match'] = idx2_original

    return cat1


input_paths = [
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_y1.fits',
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_daily.fits',
]

for input_path in input_paths:

    cat2 = Table(fitsio.read(input_path))
    print(len(cat2), len(np.unique(cat2['TARGETID'])))

    # Keep unique TARGETIDs –– choose the highest SN_B object
    cat2.sort('SN_B', reverse=True)
    _, idx_keep = np.unique(cat2['TARGETID'], return_index=True)
    cat2 = cat2[idx_keep]
    print(len(cat2), len(np.unique(cat2['TARGETID'])))

    if not np.all(cat2['REF_EPOCH']==2015.5):
        cat2['RA'] = cat2['TARGET_RA'] - (cat2['REF_EPOCH'] - 2015.5) * cat2['PMRA'] * 1e-3/3600 / np.cos(np.radians(cat2['TARGET_DEC']))
        cat2['RA'] = (cat2['RA'] + 360)%360  # Wrap around
        cat2['DEC'] = cat2['TARGET_DEC'] - (cat2['REF_EPOCH'] - 2015.5) * cat2['PMDEC'] * 1e-3/3600
    else:
        cat2['RA'], cat2['DEC'] = cat2['TARGET_RA'], cat2['TARGET_DEC']

    for field in ['north', 'south']:
        output_path = input_path.replace('.fits', '-dr9_sweep_{}.fits'.format(field)).replace('/dvs_ro/', '/global/')
        sweep_dir = os.path.join('/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field))
        cat1_paths = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))

        ra_col, dec_col, search_radius = 'RA', 'DEC', 0.1
        ra2full = np.array(cat2[ra_col])
        dec2full = np.array(cat2[dec_col])

        n_process = 256
        with Pool(processes=n_process) as pool:
            res = pool.map(do_stuff, np.arange(len(cat1_paths)))

        # Remove None elements from the list
        for index in range(len(res)-1, -1, -1):
            if res[index] is None:
                res.pop(index)

        sweep = vstack(res)
        sweep['TARGETID'] = cat2['TARGETID'][sweep['idx_match']]
        sweep.remove_column('idx_match')
        sweep['TARGETID_DR9'] = encode_targetid(sweep['OBJID'], sweep['BRICKID'], sweep['RELEASE'])

        sweep.write(output_path)

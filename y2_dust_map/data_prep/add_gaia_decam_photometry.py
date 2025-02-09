# Add Gaia-predicted LS photometry to the star catalogs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


# Coefficients for EDR3
coeffs_north = dict(
    g=[-0.1078245455, 0.3893640995, 0.5092663339, -0.3151256110,
        -0.7369068916, 1.5427860782, -0.8189261404, -0.1437262966,
        0.3244774248, -0.1454065026, 0.0316391437, -0.0034792332,
        0.0001550303],
    r=[0.1187517272, -0.2861450153, 0.0561561954, -0.0042551664,
        0.1037731685, -0.2303693495, 0.1642773825, 0.0129511209,
        -0.0667456034, 0.0338394726, -0.0079118342, 0.0009141628,
        -0.0000422475],
    z=[0.5239827217, -1.0585732272, 0.1370934028, 0.2113521602,
        -0.2514153702, 0.1056034020, -0.0193360785, 0.0013042005])
coeffs_south = dict(
    g=[-0.1189638126, 0.3175863532, 1.0168285818, -0.9809641253,
        -1.8311255772, 5.1626355784, -4.9756924458, 2.5246836016,
        -0.7372608692, 0.1225117718, -0.0101571104, 0.0002021176,
        0.0000149423],
    r=[0.1436855876, -0.2924737851, -0.0574068913, 0.1046897287,
        0.3683922665, -0.8838981047, 0.7765326441, -0.3012260862,
        0.0288217837, 0.0166817926, -0.0062362289, 0.0008462463,
        -0.0000423381],
    z=[0.5113045464, -1.0363922490, 0.1631593013, 0.1514953987,
        -0.2067508233, 0.0899108780, -0.0167211576, 0.0011369085])
bprp_min, bprp_max = -0.5, 4.5

for field in ['north', 'south']:

    cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_combined_{}.fits'.format(field),
                columns=['TARGETID', 'PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'PHOTSYS_PHOTOM']))

    # Gaia-predicted magnitudes
    for band in ['g', 'r', 'z']:
        cat['ls_mag_'+band] = np.zeros(len(cat))
        mag = np.copy(cat['PHOT_G_MEAN_MAG'])
        for photsys in ['N', 'S']:
            mask = cat['PHOTSYS_PHOTOM']==photsys
            if photsys=='N':
                coeffs = coeffs_north.copy()
            else:
                coeffs = coeffs_south.copy()
            for order, c in enumerate(coeffs[band]):
                x = (cat['PHOT_BP_MEAN_MAG']-cat['PHOT_RP_MEAN_MAG'])[mask]
                x = np.clip(x, bprp_min, bprp_max)
                mag[mask] += c * (x)**order
        cat['ls_mag_'+band] = mag

    cat = cat[['TARGETID', 'ls_mag_g', 'ls_mag_r', 'ls_mag_z']]
    cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/stars_combined_gaia_ls_photom_{}.fits'.format(field))

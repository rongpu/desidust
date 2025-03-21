# Use nside=64 instead of 128 due to the lower density of the drz stars

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
import healpy as hp

from desitarget import targetmask
import json
from sklearn.neighbors import NearestNeighbors
from multiprocessing import Pool
from scipy.stats import binned_statistic_2d

from scipy.optimize import minimize

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

field = 'south'

nmad = lambda x: 1.4826*np.median(np.abs(x-np.median(x)))
rms = lambda x: np.sqrt(np.sum(x**2)/len(x))

def rlm_fit2d(x, y, z, t=1.5, order=2):
    '''
    2D robust polynomial fit.

    Given x, y and z arrays, calculate the 2D robust
    polynomial fit of arbitrary order. Huber weight
    function is used.

    See also poly_val2d.py

    INPUT:
    1D arrays of x, y and z values; tunning parameter t;
    order of the polynomial fit.
    Note that x*y is a second order term.

    OUTPUT:
    Array of parameters of the polynomial [a0, a1, a2 ...]
    so that for an n-th order fit
    z = a0 + a1*y + a2*y**2 + ... + a_n*y**n + a_(n+1)*x
    + a_(n+2)*x*y + a_(n+3)*x*y**2 + ... + a_(2n)*x*y**(n-1)
    + a_(2n+1)*x**2 + a_(2n+2)*x**2*y + ... ...

    For example 2nd order fit gives:
    z = a0 + a1*y + a2*y**2 + a3*x + a4*x*y + a5*x**2

    New values can be evaluated using poly_val2d.py
    '''

    import statsmodels.api as sm

    ncols = (order+2)*(order+1)//2
    a = np.zeros((x.size, ncols))
    k=0
    for i in range(order+1):
        for j in range(order-i+1):
            a[:,k] = x**i * y**j
            k+=1
    res = sm.RLM(z, a, M=sm.robust.norms.HuberT(t=t)).fit()
    m = res.params
    return(m)


def poly_val2d(x, y, m):
    '''
    Evaluate the 2D polynomial from x and yvalues and polynomial
    parameters

    See also rlm_fit2d.py

    INPUT:
    1D array of x and y values;
    1D array of polynomial parameters (generated by rlm_fit2d.py).

    OUTPUT:
    1D array of the evaluated values of the polynomial.
    '''

    order = int((np.sqrt(8*len(m)+1)-3)/2)
    z = np.zeros(x.shape)
    k=0
    for i in range(order+1):
        for j in range(order-i+1):
            z += m[k] * x**i * y**j
            k+=1
    return z

cat = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/rz_corrected/rz_{}.fits'.format(field)))
print(len(cat))

cat_all = cat.copy()

nside = 64
cat['HPXPIXEL'] = hp.ang2pix(nside, cat['TARGET_RA'], cat['TARGET_DEC'], nest=False, lonlat=True)

maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/maps/drz_map_{}_{}.fits'.format(field, nside)))

mask = np.in1d(maps['HPXPIXEL'], cat['HPXPIXEL'])
print(np.sum(mask)/len(mask))
maps = maps[mask]

mask = maps['n_star']>=32
print('n_star', np.sum(mask)/len(mask))
maps = maps[mask]

cat = join(cat, maps[['delta_rz_mean', 'delta_rz_hlmean', 'HPXPIXEL', 'n_star']], keys='HPXPIXEL', join_type='inner')

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

mask = cat['EBV_SFD']<0.1
print('EBV_SFD', np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]
print(len(cat))

cat['drz_diff_raw'] = cat['drz_raw'] - cat['delta_rz_mean']
cat['drz_diff'] = cat['drz'] - cat['delta_rz_mean']

min_rz_err = 0.012
cat_all['rz_err'] = 0.

# star_type = 'bright_mws'
for star_type in ['dark_std', 'bright_std', 'bright_mws', 'backup_std', 'backup_mws']:

    print(star_type)
    mask0 = cat[star_type].copy()
    print(np.sum(mask0), np.sum(mask0)/len(mask0))

    # compute the average r-z error in 2D bins of SN_B and model_rz_std; count the number of objects in each bin
    xbins, ybins = np.linspace(0, 80, 81), np.linspace(0.01, 0.036, 76)
    ret = binned_statistic_2d(cat['SN_B'][mask0], cat['model_rz_std'][mask0], cat['drz_diff'][mask0], statistic=rms, bins=[xbins, ybins])
    count = binned_statistic_2d(cat['SN_B'][mask0], cat['model_rz_std'][mask0], None, statistic='count', bins=[xbins, ybins])

    # get the binned_statistic_2d bin indices
    # (for determining which regions in the SN_B and model_rz_std are removed by the minimum count cut)
    xx, yy = np.meshgrid((xbins[1:]+xbins[:-1])/2, (ybins[1:]+ybins[:-1])/2)
    tmp = binned_statistic_2d(xx.flatten(), yy.flatten(), None, statistic='count', bins=[xbins, ybins])
    bin_idx = tmp.binnumber

    # require a minimum count
    mask_count = count.statistic>=16
    ret.statistic[~mask_count] = np.nan
    bin_idx = bin_idx[mask_count.T.flatten()]

    rms_data = ret.statistic.T.flatten()
    sn = xx.flatten()
    model_std = yy.flatten()

    # polynomial fit
    mask = ~np.isnan(rms_data)
    coeffs = rlm_fit2d(sn[mask], model_std[mask], rms_data[mask], t=0.01, order=5)
    # print(coeffs)

    # counts for the full catalog
    count1 = binned_statistic_2d(cat['SN_B'], cat['model_rz_std'], cat['drz_diff'], statistic='count', bins=[xbins, ybins])
    mask1 = mask0 & np.in1d(count1.binnumber, bin_idx)
    print('Outliers in SN_B vs model_rz_std: {:.2f}%'.format(100*np.sum(~mask1 & mask0)/np.sum(mask0)))

    # compute the error for the well-sampled regions in SN_B vs model_rz_std space
    count_all = binned_statistic_2d(cat_all['SN_B'], cat_all['model_rz_std'], None, statistic='count', bins=[xbins, ybins])
    good_region = cat_all[star_type] & np.in1d(count_all.binnumber, bin_idx)
    cat_all['rz_err'][good_region] = np.clip(poly_val2d(cat_all['SN_B'][good_region], cat_all['model_rz_std'][good_region], coeffs), min_rz_err, np.inf)

    # create nearest-neighbor training data
    x1, x2 = cat_all['SN_B'].copy(), cat_all['model_rz_std'].copy()
    x1_sigma, x2_sigma = nmad(x1), nmad(x2)
    x1 /= x1_sigma
    x2 /= x2_sigma
    X = np.vstack([x1[good_region], x2[good_region]]).T
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(X)

    # Compute the error for the "outlier" regions
    mask = cat_all[star_type] & (~np.in1d(count_all.binnumber, bin_idx))
    X_outliers = np.vstack([cat_all['SN_B'][mask]/x1_sigma, cat_all['model_rz_std'][mask]/x2_sigma]).T
    idx = nbrs.kneighbors(X_outliers)[1]
    cat_all['rz_err'][mask] = cat_all['rz_err'][good_region][idx[:, 0]]
    assert np.sum(cat_all['rz_err'][cat_all[star_type]]==0)==0

    ############################################# Make plots #############################################

    rms_predict = poly_val2d(sn, model_std, coeffs)
    mask = ~np.isnan(rms_data)
    rms_predict[~mask] = np.nan

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    pcm = ax.imshow(count.statistic.T, origin='lower', extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], aspect='auto', vmin=0, vmax=np.percentile(count.statistic, 99.5), cmap='gray_r')
    ax.set_xlabel('SN_B')
    ax.set_ylabel('model_rz_std')
    fig.colorbar(pcm, ax=ax)
    ax = axes[1]
    pcm = ax.imshow(rms_data.reshape(ret.statistic.T.shape), origin='lower', extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], aspect='auto', vmin=0.01, vmax=0.035, cmap='jet')
    ax.set_xlabel('SN_B')
    ax.set_ylabel('model_rz_std')
    fig.colorbar(pcm, ax=ax)
    plt.show()

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    pcm = ax.imshow((rms_predict-rms_data).reshape(ret.statistic.T.shape), origin='lower', extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], aspect='auto', vmin=-0.01, vmax=0.01, cmap='seismic')
    ax.set_xlabel('SN_B')
    ax.set_ylabel('model_rz_std')
    fig.colorbar(pcm, ax=ax)
    ax = axes[1]
    pcm = ax.imshow(rms_predict.reshape(ret.statistic.T.shape), origin='lower', extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], aspect='auto', vmin=0.01, vmax=0.035, cmap='jet')
    ax.set_xlabel('SN_B')
    ax.set_ylabel('model_rz_std')
    fig.colorbar(pcm, ax=ax)
    plt.show()

    count2 = binned_statistic_2d(cat['SN_B'][mask1], cat['model_rz_std'][mask1], None, statistic='count', bins=[xbins, ybins])
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    pcm = ax.imshow(count2.statistic.T, origin='lower', extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], aspect='auto', vmin=0, vmax=1, cmap='gray_r')
    ax.set_xlabel('SN_B')
    ax.set_ylabel('model_rz_std')
    fig.colorbar(pcm, ax=ax)
    ax = axes[1]
    mask = cat_all[star_type].copy()
    ret1 = binned_statistic_2d(cat_all['SN_B'][mask], cat_all['model_rz_std'][mask], cat_all['rz_err'][mask], statistic=np.mean, bins=[xbins, ybins])
    pcm = ax.imshow(ret1.statistic.T, origin='lower', extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()], aspect='auto', vmin=0.01, vmax=0.035, cmap='jet')
    ax.set_xlabel('SN_B')
    ax.set_ylabel('model_rz_std')
    fig.colorbar(pcm, ax=ax)
    plt.show()

    plt.hist(rms_data, 100, range=(0.01, np.nanpercentile(rms_data, 99.)+0.005), alpha=0.5)
    plt.hist(rms_predict, 100, range=(0.01, np.nanpercentile(rms_data, 99.)+0.005), alpha=0.5)
    plt.hist(rms_predict-rms_data + np.nanmean(rms_data), 100,
    #              range=(np.nanmin(rms_data), np.nanmax(rms_data)), alpha=0.5, histtype='step', color='r', lw=1)
             range=(np.nanmin(rms_data), 0.035), alpha=0.5, histtype='step', color='r', lw=1)
    plt.show()

assert np.sum(cat_all['rz_err']==0)==0

print(len(cat_all), len(np.unique(cat_all['TARGETID_DR9'])), len(cat_all)==len(np.unique(cat_all['TARGETID_DR9'])))

cat_all = cat_all[['TARGETID_DR9', 'rz_err']]
cat_all.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/rz_corrected/rz_{}_predicted_error.fits'.format(field), overwrite=True)




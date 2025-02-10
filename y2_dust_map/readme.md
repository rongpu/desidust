Code for the DESI dust map https://arxiv.org/abs/2409.05140

## Data prep
 - Create combined coadd and model spectrum files (for reproducibility and archiving)
    `data_prep/combine_coadd_and_specmod_spectra-y1.py`
    `data_prep/combine_coadd_and_specmod_spectra-daily.py`
 - Create catalogs with synthetic photometry
    `data_prep/create_star_catalog-y1.py`
    `data_prep/create_star_catalog-daily.py`
 - Match to LS DR9 and Gaia DR3
    `data_prep/match_to_dr9_sweeps.py`
    `data_prep/match_to_gaia_dr3.py`
 - Combine the catalogs
    `data_prep/combine_catalogs.py`
 - Add Gaia-predicted DECam photometry
    `data_prep/add_gaia_decam_photometry.py`

## Create delta_g-r catalogs and maps
 - Create reference catalog for g-r correction
    `create_dgr_reference_catalog.py`
 - Create catalogs and unweighted maps with corrected delta_g-r
    `create_dgr_catalog_and_unweighted_maps.py`
 - Compute the per-star errors for optimal weighting in map creation
    `per_star_weights_dgr.py`
 - Create separate North and South weighted maps (e.g., for determining the North/South relative zero point offset)
    `create_ns_dgr_weighted_maps.py`
 - Created final North+South combined catalog and maps; the North/South relative zero point offset is corrected for; extra cuts (parallax and SN_B) are applied on the final catalog
    `create_combined_dgr_catalogs_and_weighted_maps.py`

## Create delta_r-z catalogs and maps
 - Create reference catalog for r-z correction
    `create_drz_reference_catalog.py`
 - Create catalogs and unweighted maps with corrected delta_r-z
    `create_drz_catalog_and_unweighted_maps.py`
 - Compute the per-star errors for optimal weighting in map creation
    `per_star_weights_drz.py`
 - Create separate North and South weighted maps (e.g., for determining the North/South relative zero point offset)
    `create_ns_drz_weighted_maps.py`
 - Created final North+South combined catalog and maps; the North/South relative zero point offset is corrected for; extra cuts (parallax and SN_B) are applied on the final catalog
    `create_combined_drz_catalogs_and_weighted_maps.py`

## Create the KP3 map and the final public map
 - Create smoothed and inpainted maps
    `create_smoothed_maps.py`
 - Create KP3 maps; it uses DECam extinction coefficients for both North and South
    `create_v1_dust_map_for_kp3.py`
 - Create "fill frac" maps
    `create_fill_frac_maps.py`
 - Create final public maps released with the paper; it uses different sets of extinction coefficients for North and South
    `create_final_maps.py`

## Misc
 - Create smoothed and inpainted maps for SFD
    `create_smoothed_sfd_maps.py`


## Dependencies
Most of the code is run with:
```
numpy 1.22.4
matplotlib 3.6.2
astropy 5.2.1
fitsio 1.1.8
sklearn 1.1.1
desitarget 2.6.0
scipy 1.8.1
statsmodels 0.13.5
speclite 0.16
desispec 0.56.5
desiutil 3.2.6
```
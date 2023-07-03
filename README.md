# Hydrodynamical effects on higher-order lensing statistics

This code works out the impact of hydrodynamical baryonic feedback effects on 2-, 3 and 4-point weak lensing data statistics. This was used in the numerical analysis of the paper:

- Barreira et al 2020, https://arxiv.org/abs/1904.02070

The calculation of the 3- and 4-point functions uses the formalism of the (Response Approach to Perturbation Theory)[https://arxiv.org/abs/1703.09212].

### Dependencies

- python (numpy, scipy, matplotlib)

### Code overview

- prepare_for_lensing.py
  - This defines global parameters and functions. It is imported by the other scripts.

- compute_kappa_spectra.py
  - This computes the 2-point and the squeezed 3-point function
  
- compute_cov_g.py and compute_cov_ssc_limber.py
  - This computes the 2-point function covariance (which is a 4-point function). It includes the two most important contributions: the Gaussian (cov_g) and the super-sample covariance (cov_ssc) terms.

- lookup_tables/
  - This containts tables with 3D power spectra, as well as N-body simulation response function data used in the calculation

- plots/
  - The plotting scripts are in this folder

### Gallery

Summary of the kinematic regimes in $k_1-k_2$ space

<img src="plots/.png" width="1000" height=auto/>


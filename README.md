# Hydrodynamical effects on higher-order lensing statistics

This code works out the impact of hydrodynamical baryonic feedback effects on 2-, 3 and 4-point weak lensing data statistics. This was used in the numerical analysis of the paper:

- Barreira et al 2020, https://arxiv.org/abs/1904.02070

The calculation of the 3- and 4-point functions uses the formalism of the (Response Approach to Perturbation Theory)[https://arxiv.org/abs/1703.09212].

### Dependencies

- python (numpy, scipy, matplotlib)

### Code overview

- The code prepare_for_lensing.py defines global parameters and functions. It is imported by the other scripts.

- The file compute_kappa_spectra.py computes the 2-point and the squeezed 3-point function.
  
- The files compute_cov_g.py and compute_cov_ssc_limber.py compute the 2-point function covariance (which is a 4-point function). It includes the two most important contributions: the Gaussian (cov_g) and the super-sample covariance (cov_ssc) terms.

- The folder lookup_tables/ contains tables with 3D power spectra, as well as N-body simulation response function data used in the calculation.

- The folder plots/ contains the plotting scripts.

### Gallery

Summary of the kinematic regimes in $k_1-k_2$ space

<img src="plots/.png" width="1000" height=auto/>


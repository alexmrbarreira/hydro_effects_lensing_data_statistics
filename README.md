# Hydrodynamical effects on higher-order lensing statistics

This code works out the impact of hydrodynamical baryonic feedback effects on 2-, 3 and 4-point weak lensing data statistics. This was used in the numerical analysis of the paper:

- Barreira et al 2020, https://arxiv.org/abs/1904.02070

The calculation of the 3- and 4-point functions uses the formalism of the (Response Approach to Perturbation Theory)[https://arxiv.org/abs/1703.09212].

### Dependencies

- python (numpy, scipy, matplotlib)

### Code overview

- prepare_for_lensing.py

This 

- The scripts in compute_cov/ execute the covariance calculation
- The scripts in plots/ make plots (figures are stored here too).

### Gallery

Summary of the kinematic regimes in $k_1-k_2$ space

<img src="plots/fig_regimes_v2.png" width="1000" height=auto/>


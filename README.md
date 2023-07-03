# Hydrodynamical effects on higher-order lensing statistics

This code works out the impact of hydrodynamical baryonic feedback effects on 2-, 3 and 4-point weak lensing statistics. This was used in the numerical analysis of the paper:

- Barreira et al 2020, [https://arxiv.org/abs/1904.02070]

who developed the response approach and used it, as a first application, to address the problem of accurate analytical calculations of covariance matrices. This code includes the Gaussian (G) term and the connected non-Gaussian (cNG) term up to 1-loop. For the super-sample covariance (SSC) term see [this other code](https://github.com/alexmrbarreira/gravitational_lensing_SSC).

For more applications of the response approach see:

- Barreira, Krause & Schmidt 2018, https://arxiv.org/abs/1711.07467
- Barreira 2019, https://arxiv.org/abs/1901.01243
- Barreira et al. 2019, https://arxiv.org/abs/1904.02070
- Halder & Barreira 2022, https://arxiv.org/abs/2201.05607

### Dependencies

- python (numpy, scipy, matplotlib)

### Code overview

- The files parameters.py and functions.py define global parameters, variables and functions
- The scripts in compute_cov/ execute the covariance calculation
- The scripts in plots/ make plots (figures are stored here too).

### Gallery

Summary of the kinematic regimes in $k_1-k_2$ space

<img src="plots/fig_regimes_v2.png" width="1000" height=auto/>

Stitching of the tree-level standard perturbation theory (SPT) and response results

<img src="plots/fig_tree_transition_mono.png" width="1000" height=auto/>

Contributions along the diagonal

<img src="plots/fig_covariance_diagonal_mono.png" width="600" height=auto/>

Response functions used in the calculation

<img src="plots/fig_responses_eul.png" width="1000" height=auto/>

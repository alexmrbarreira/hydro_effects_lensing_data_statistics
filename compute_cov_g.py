from prepare_for_lenscov import *

# ================================================================================ #
# Load power spectra 
# ================================================================================ #

C_l_dmo = load('data_store/C_l_dmo.npy')
C_l_tng = load('data_store/C_l_tng.npy')

# ================================================================================ #
# Compute Gaussian covariance
# ================================================================================ #
print ('Computing Gaussian covariance ... ')

Clnoise  = sigma_eps**2./sounumden/2.

cov_l1l2_dmo = zeros([len(ll_hard), len(ll_hard)])
cov_l1l2_tng = zeros([len(ll_hard), len(ll_hard)])

for i in range(len(ll_hard)):
    deltal_now           = deltal_array[i]
    l_now                = ll_hard[i]

    # DMO
    cov_l1l2_dmo[i,i] = (4.*pi/Omega_S/l_now/deltal_now) * (C_l_dmo[i] + Clnoise)**2.
    # Hydro
    cov_l1l2_tng[i,i] = (4.*pi/Omega_S/l_now/deltal_now) * (C_l_tng[i] + Clnoise)**2.


# ================================================================================ #
# Compute Gaussian covariance
# ================================================================================ #
print ('Saving ... ')

save('data_store/cov_g_dmo', cov_l1l2_dmo)
save('data_store/cov_g_tng', cov_l1l2_tng)


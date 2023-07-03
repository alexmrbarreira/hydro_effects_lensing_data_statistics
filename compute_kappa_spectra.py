from prepare_for_lenscov import *

# ================================================================================ #
# Compute power spectrum 
# ================================================================================ #
print ('Computing kappa power spectrum ... ')
C_l_dmo = zeros(len(ll_hard))
C_l_tng = zeros(len(ll_hard))

for m in range(len(ll_hard)):
    k1_now     = (ll_hard[m])/chi_array

    # Gravity-only calculation
    Pnl_dmo           = diagonal(flipud(array(Pnl_int(z_array, k1_now))))
    C_l_dmo_integrand = Pnl_dmo * gkernel_array**2. / chi_array**2.
    C_l_dmo[m]        = integrate.trapz(C_l_dmo_integrand, chi_array)

    # Hydro calculation
    Pnl_tng           = diagonal(flipud(array(Pnl_tng_int(z_array, k1_now))))
    C_l_tng_integrand = Pnl_tng * gkernel_array**2. / chi_array**2.
    C_l_tng[m]        = integrate.trapz(C_l_tng_integrand, chi_array)

# ================================================================================ #
# Compute bispectrum
# ================================================================================ #
print ('Computing kappa bispectrum ... ')
B_lll_dmo = zeros([len(ll_hard), len(ll_soft)])
B_lll_tng = zeros([len(ll_hard), len(ll_soft)])

for i in range(len(ll_soft)):
    k_soft_now = ll_soft[i]/chi_array
    for j in range(len(ll_hard)):
        k_hard_now = ll_hard[j]/chi_array

        # Gravity only calculation
        P_hard_now_dmo      = diagonal(flipud(array(Pnl_int(z_array, k_hard_now)))) 
        P_soft_now_dmo      = diagonal(flipud(array(Plin_int(z_array, k_soft_now)))) 
        R1_now_dmo          = diagonal(flipud(array(R_1_tng_dmo_int(z_array, k_hard_now))))
        RK_now_dmo          = diagonal(flipud(array(R_K_int(z_array, k_hard_now))))
        B_lll_dmo_integrand = (R1_now_dmo + RK_now_dmo/6.) * P_hard_now_dmo * P_soft_now_dmo * gkernel_array**3. / chi_array**4. #1/6 comes from 2D angle average of [mu^2 - 1/3]
        B_lll_dmo[j,i]      =  integrate.trapz(B_lll_dmo_integrand, chi_array)

        # Hydro calculation
        P_hard_now_tng      = diagonal(flipud(array(Pnl_tng_int(z_array, k_hard_now))))
        P_soft_now_tng      = diagonal(flipud(array(Plin_int(z_array, k_soft_now))))
        R1_now_tng          = diagonal(flipud(array(R_1_tng_int(z_array, k_hard_now))))
        RK_now_tng          = diagonal(flipud(array(R_K_tng_int(z_array, k_hard_now))))
        B_lll_tng_integrand = (R1_now_tng + RK_now_tng/6.) * P_hard_now_tng * P_soft_now_tng * gkernel_array**3. / chi_array**4. #1/6 comes from 2D angle average of [mu^2 - 1/3]
        B_lll_tng[j,i]      =  integrate.trapz(B_lll_tng_integrand, chi_array)


print ('Saving ...')

save('data_store/ll_hard', ll_hard)
save('data_store/ll_soft', ll_soft)
save('data_store/C_l_dmo', C_l_dmo)
save('data_store/C_l_tng', C_l_tng)
save('data_store/B_lll_dmo', B_lll_dmo)
save('data_store/B_lll_tng', B_lll_tng)


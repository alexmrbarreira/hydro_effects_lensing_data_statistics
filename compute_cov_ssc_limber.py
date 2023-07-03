from prepare_for_lenscov import *

# ================================================================================ #
# Do variance integrals 
# ================================================================================ #
print ('Doing variance integrals ... ; same for Hydro and Gravity-only')

ll_window = 10.**linspace(-6., log10(790.), 10000)

varint_chi = zeros(len(chi_array))

for i in range(len(chi_array)):
    theta_S       = sqrt(Omega_S/pi)
    z_now         = zofchi_int(chi_array[i])
    integrand_now = ll_window**2. * (2.*special.j1(ll_window*theta_S)/ll_window/theta_S)**2. * Plin_int(z_now, ll_window/chi_array[i])[:,0]
    
    varint_chi[i] = integrate.trapz(integrand_now, log(ll_window))/2./pi

# ================================================================================ #
# Do chi integrals 
# ================================================================================ #
print ('Doing chi integrals ... ')

cov_l1l2_dmo = zeros([len(ll_hard), len(ll_hard)])
cov_l1l2_tng = zeros([len(ll_hard), len(ll_hard)])

# Loop over l1 and l2 values
for m in range(len(ll_hard)):
    print (m, 'of', len(ll_hard), '; it gets faster though ... ')
    for n in range(m, len(ll_hard)):
        k1_now = (ll_hard[m])/chi_array
        k2_now = (ll_hard[n])/chi_array

        # DMO
        R_1_l1 = diagonal(flipud(array(R_1_tng_dmo_int(z_array, k1_now))))
        R_1_l2 = diagonal(flipud(array(R_1_tng_dmo_int(z_array, k2_now))))
        R_K_l1 = diagonal(flipud(array(R_K_int(z_array, k1_now))))
        R_K_l2 = diagonal(flipud(array(R_K_int(z_array, k2_now))))
        Pnl_l1 = diagonal(flipud(array(Pnl_int(z_array, k1_now))))
        Pnl_l2 = diagonal(flipud(array(Pnl_int(z_array, k2_now))))
        ssc_integrand     =  varint_chi * (R_1_l1 + R_K_l1/6.) * (R_1_l2 + R_K_l2/6.) * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**6.
        cov_l1l2_dmo[m,n] = integrate.trapz(ssc_integrand, chi_array)

        # Hydro 
        R_1_l1 = diagonal(flipud(array(R_1_tng_int(z_array, k1_now))))
        R_1_l2 = diagonal(flipud(array(R_1_tng_int(z_array, k2_now))))
        R_K_l1 = diagonal(flipud(array(R_K_tng_int(z_array, k1_now))))
        R_K_l2 = diagonal(flipud(array(R_K_tng_int(z_array, k2_now))))
        Pnl_l1 = diagonal(flipud(array(Pnl_tng_int(z_array, k1_now))))
        Pnl_l2 = diagonal(flipud(array(Pnl_tng_int(z_array, k2_now))))
        ssc_integrand     =  varint_chi * (R_1_l1 + R_K_l1/6.) * (R_1_l2 + R_K_l2/6.) * Pnl_l1 * Pnl_l2 * gkernel_array**4. / chi_array**6.
        cov_l1l2_tng[m,n] = integrate.trapz(ssc_integrand, chi_array)

# ================================================================================ #
# Symmetrize and save 
# ================================================================================ #
print ('Symmetrizing and saving ... ')

cov_l1l2_dmo = symmetrize_matrix(cov_l1l2_dmo)
cov_l1l2_tng = symmetrize_matrix(cov_l1l2_tng)

save('data_store/cov_ssc_limber_dmo', cov_l1l2_dmo)
save('data_store/cov_ssc_limber_tng', cov_l1l2_tng)



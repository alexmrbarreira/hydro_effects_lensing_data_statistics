import sys; sys.path.append('../'); 
from prepare_for_lenscov import *

## ========================================================
## Get what to plot
## ========================================================

data_vector_tng = load('../data_store/C_l_tng.npy')
data_vector_dmo = load('../data_store/C_l_dmo.npy')
matrix_C1C2_tng = outer(data_vector_tng, data_vector_tng)
matrix_C1C2_dmo = outer(data_vector_dmo, data_vector_dmo)

# Load Limber SSC
cov_ssc_dmo = load('../data_store/cov_ssc_limber_dmo.npy')
cov_ssc_tng = load('../data_store/cov_ssc_limber_tng.npy')

# Load G
cov_g_dmo = load('../data_store/cov_g_dmo.npy')
cov_g_tng = load('../data_store/cov_g_tng.npy')

# Make a total covariance and compute their inverses
cov_dmo = cov_g_dmo + cov_ssc_dmo
cov_tng = cov_g_tng + cov_ssc_tng

invcov_dmo = linalg.inv(cov_dmo)
invcov_tng = linalg.inv(cov_tng)

## ========================================================
## Compute cumulative Signal-to-Noise
## ========================================================
# Function that computes the cumulative Signal-to-Noise ratio as a function of k_max
def compute_cumSN_powerspectrum(x, y, cov):
    # Arrays that'll store kmax values and corresponding S/N
    xmax_array = zeros(len(x))
    sn_array   = zeros(len(x))
    indices_to_delete = [] # will store indices to tell which columns/lines/entries to progressively remove
    # Loop over number of k_bins
    for i in range(len(x)):
        # Delete wanted number of entries
        cov_now = delete(cov    , indices_to_delete, 1)
        cov_now = delete(cov_now, indices_to_delete, 0)
        x_now   = delete(x      , indices_to_delete)
        y_now   = delete(y      , indices_to_delete)
        # Invert the resulting reduced covariance matrix
        invcov_now     = linalg.inv(cov_now)
        # Compute the S/N (note the array is filled backwards)
        xmax_array[len(x) - 1 - i] = max(x_now)
        sn_array[len(x) - 1 - i]   = sqrt(dot(y_now, dot(invcov_now  , y_now)))
        # Increment number of indices to delete
        indices_to_delete  = [len(x) - 1 - j for j in range(i+1)]
    # Return
    return xmax_array, sn_array

# Note that ellmax is the same for all
ellmax, snratio_g_dmo       = compute_cumSN_powerspectrum(ll_hard, data_vector_tng, cov_g_dmo) 
ellmax, snratio_g_tng       = compute_cumSN_powerspectrum(ll_hard, data_vector_tng, cov_g_tng) 
ellmax, snratio_g_ssc_dmo   = compute_cumSN_powerspectrum(ll_hard, data_vector_tng, cov_g_dmo + cov_ssc_dmo)
ellmax, snratio_g_ssc_tng   = compute_cumSN_powerspectrum(ll_hard, data_vector_tng, cov_g_tng + cov_ssc_tng) 

## ====================================================================================================
## Draw samples with TNG mean and TNG covariance and compute chi^2 for varying DMO/TNG models and covs.
## ====================================================================================================

Nrandom = 1000

chi2_tng_tng = zeros(Nrandom) #chi^2 for TNG model and TNG covariance
chi2_tng_dmo = zeros(Nrandom) #chi^2 for TNG model and DMO covariance
chi2_dmo_tng = zeros(Nrandom) #chi^2 for DMO model and TNG covariance

for i in range(Nrandom):
    sample_data_now = random.multivariate_normal(data_vector_tng, cov_tng)

    chi2_tng_tng[i] = dot( dot( (data_vector_tng - sample_data_now), invcov_tng ) , (data_vector_tng - sample_data_now) )
    chi2_tng_dmo[i] = dot( dot( (data_vector_tng - sample_data_now), invcov_dmo ) , (data_vector_tng - sample_data_now) )
    chi2_dmo_tng[i] = dot( dot( (data_vector_dmo - sample_data_now), invcov_tng ) , (data_vector_dmo - sample_data_now) )

## ========================================================
## Plot
## ========================================================
labelsize = 26
ticksize  = 26
textsize  = 30
titlesize = 24
text_font = 22
legend_font = 22
v_min = 0.01
v_max = 0.60

fig3 = plt.figure(3, figsize=(17., 7.))
fig3.subplots_adjust(left=0.055, right=0.985, top=0.93, bottom=0.14, wspace = 0.25, hspace = 0.35)

# Plot covariance ratio hydro/gravity
panel = fig3.add_subplot(1,2,1)
# Add the data
pcolor(ll_hard, ll_hard, (cov_g_tng + cov_ssc_tng)/(cov_g_dmo + cov_ssc_dmo), vmin = 0.80, vmax = 1., cmap = 'jet')
# Cosmetics
panel.set_aspect('equal')
colorbar().ax.tick_params(labelsize=ticksize)
annotate(r"Hydro/Gravity", xy = (0.04, 0.90), xycoords='axes fraction',
xytext = (0.04, 0.90), textcoords='axes fraction', bbox=dict(boxstyle="round4", fc = 'w', alpha = 0.6), color = 'k', fontsize = textsize-6)
title(r'${\rm Cov}^{\rm G} + {\rm Cov}^{\rm SSC}$'   , fontsize = titlesize)
xlabel(r'$\ell_1$' , fontsize = labelsize)
ylabel(r'$\ell_2$' , fontsize = labelsize)
xscale('log')
yscale('log')
xlim(min(ll_hard), max(ll_hard))
ylim(min(ll_hard), max(ll_hard))
xticks(size = ticksize)
yticks(size = ticksize)

# Plot the signal-to-noise ratio hydro/gravity
fig3.add_subplot(2,2,2)
# Add the data
plot(ellmax, snratio_g_ssc_tng/snratio_g_ssc_dmo, linewidth = 1.5, linestyle = 'solid' , c = 'r', label = r'Hydro/Gravity')
# Cosmetics
axhline(1, linewidth = 1.5, linestyle = 'dotted', c = 'k')
xlabel(r'$\ell_{\rm max}$' , fontsize = labelsize+2)
ylabel(r'S/N$(<\ell_{\rm max})$'        , fontsize = labelsize+2)
xscale('log')
xlim(min(ll_hard), max(ll_hard))
ylim(0.99, 1.06)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font-1.5}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 1)

# Plot impact on the goodness of fit from ignoring hydro effects
fig3.add_subplot(2,2,4)
# Add the data
toplot_1 =     (chi2_tng_dmo - chi2_tng_tng)/len(ll_hard)
toplot_2 =     (chi2_dmo_tng - chi2_tng_tng)/len(ll_hard)
valuemin = 5.
toplot_1_new = valuemin * toplot_1/abs(min(toplot_1))
hist(toplot_1_new, bins = 50, linewidth = 2., color = 'k', histtype='step', label = r'No\ baryons\ in\ covariance')
hist(toplot_2    , bins = 50, linewidth = 2., color = 'r', histtype='step', label = r'No\ baryons\ in\ signal')
# Cosmetics
xlabel(r'$\Delta \chi^2 / dof$' , fontsize = labelsize+2)
ylabel(r'Counts'        , fontsize = labelsize+2)
nticks_1 = 2
fake_tick_values = linspace(min(toplot_1_new), max(toplot_1_new), nticks_1).tolist() + [5, 10, 15, 20, 25, 30]
real_tick_values = []
real_tick_number = linspace(min(toplot_1)    , max(toplot_1)    , nticks_1)
for i in range(len(real_tick_number)):
    real_tick_values.append(str("%.2f" % real_tick_number[i]))
real_tick_values  = real_tick_values + ["5", "10", "15", "20", "25", "30"]
axvline(max(toplot_1_new), linewidth = 1, linestyle='dashed', c='k')
xticks(fake_tick_values, real_tick_values, size = ticksize)
yticks(size = ticksize)
xlim(-valuemin-1., 27)
params = {'legend.fontsize': legend_font-4}; pylab.rcParams.update(params); legend(loc = 'upper center', ncol = 1)

fig3.savefig('fig_covariance.png')

#show()

import sys; sys.path.append('../'); 
from prepare_for_lenscov import *

## ========================================================
## Load arrays to plot
## ========================================================

C_l_dmo   = load('../data_store/C_l_dmo.npy')
C_l_tng   = load('../data_store/C_l_tng.npy')
B_lll_dmo = load('../data_store/B_lll_dmo.npy')
B_lll_tng = load('../data_store/B_lll_tng.npy')

cov_ssc_dmo = load('../data_store/cov_ssc_limber_dmo.npy')
cov_ssc_tng = load('../data_store/cov_ssc_limber_tng.npy')

## ========================================================
## Plot them 
## ========================================================

labelsize = 26
ticksize  = 26
textsize  = 30
titlesize = 26
text_font = 18
legend_font = 22

lfac = ll_hard*(ll_hard+1.) / 2. / pi
fsq      = 5.

fig4 = plt.figure(4, figsize=(8., 10.))
fig4.subplots_adjust(left=0.165, right=0.98, top=0.95, bottom=-0.20, hspace = 0.15, wspace = 0.30)

# Plot the spectra
fig4.add_subplot(2,1,1)
# Add the data
plot(ll_hard, lfac       *C_l_dmo*1.0e30, linewidth = 1.5, linestyle = 'dashed', c = 'k', label = r'Gravity')
plot(ll_hard, lfac       *C_l_tng*1.0e30, linewidth = 1.5, linestyle = 'solid' , c = 'k', label = r'Hydro')
plot(ll_hard, lfac       *C_l_dmo, linewidth = 1.5, linestyle = 'dashed', c = 'b')
plot(ll_hard, lfac       *C_l_tng, linewidth = 1.5, linestyle = 'solid' , c = 'b')
wherecond = where(ll_hard >= fsq*ll_soft[0])
plot(ll_hard[wherecond], (1.0e5*lfac*B_lll_dmo[:,0])[wherecond], linewidth = 1.5, linestyle = 'dashed', c = 'g')
plot(ll_hard[wherecond], (1.0e5*lfac*B_lll_tng[:,0])[wherecond], linewidth = 1.5, linestyle = 'solid' , c = 'g')
wherecond = where(ll_hard >= fsq*ll_soft[2])
plot(ll_hard[wherecond], (1.0e5*lfac*B_lll_dmo[:,2])[wherecond], linewidth = 1.5, linestyle = 'dashed', c = 'r')
plot(ll_hard[wherecond], (1.0e5*lfac*B_lll_tng[:,2])[wherecond], linewidth = 1.5, linestyle = 'solid' , c = 'r')
plot(ll_hard, lfac*1.0e-3*diagonal(cov_ssc_dmo)/C_l_tng**2., linewidth = 1.5, linestyle = 'dashed', c = 'm')
plot(ll_hard, lfac*1.0e-3*diagonal(cov_ssc_tng)/C_l_tng**2., linewidth = 1.5, linestyle = 'solid' , c = 'm')
# Cosmetics
annotate(r"$10^{-3} \times {\rm Cov}^{SSC}_\kappa(\ell_h, \ell_h) / C_\kappa(\ell_h)^2$"  , xy = (0.30,0.46), xycoords='axes fraction', color = 'm', fontsize = text_font, rotation = 40.)
annotate(r"$C_\kappa(\ell_h)$"                                                            , xy = (0.50,0.35), xycoords='axes fraction', color = 'b', fontsize = text_font, rotation = 18.)
annotate(r"$10^5 \times \mathcal{B}_\kappa(\ell_h, \ell_h, \ell_s = 20)$"                 , xy = (0.50,0.21), xycoords='axes fraction', color = 'g', fontsize = text_font, rotation = 18.)
annotate(r"$10^5 \times \mathcal{B}_\kappa(\ell_h, \ell_h, \ell_s = 100)$"                , xy = (0.50,0.05), xycoords='axes fraction', color = 'r', fontsize = text_font, rotation = 18.)
annotate(r"$z_s = 1$"                                                                     , xy = (0.03,0.65), xycoords='axes fraction', color = 'k', fontsize = text_font+8)
ylabel(r'$\{C_\kappa, \mathcal{B}_\kappa, {\rm Cov}^{SSC}_\kappa\}\times\ \ell\left(\ell + 1\right) / (2\pi)$'        , fontsize = labelsize+2)
xscale('log')
yscale('log')
xlim(min(ll_hard), max(ll_hard))
ylim(1.0e-6, 3.0e-2)
xticks(size = ticksize+1)
yticks(size = ticksize+1)
params = {'legend.fontsize': legend_font-1.5}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 1)

# Plot ratio of hydro/gravity-only result
fig4.add_subplot(4,1,3)
# Add the data
plot(ll_hard,   C_l_tng     /  C_l_dmo                               , linewidth = 1.5, linestyle = 'solid', c = 'b', label = r'$C_\kappa(\ell_h)$')
wherecond = where(ll_hard >= fsq*ll_soft[0])
plot(ll_hard[wherecond], (B_lll_tng[:,0]/B_lll_dmo[:,0])[wherecond]  , linewidth = 1.5, linestyle = 'solid', c = 'g', label = r'$\mathcal{B}_\kappa(\ell_h, \ell_h, \ell_s = 20)$')
wherecond = where(ll_hard >= fsq*ll_soft[2])
plot(ll_hard[wherecond], (B_lll_tng[:,2]/B_lll_dmo[:,2])[wherecond]  , linewidth = 1.5, linestyle = 'solid', c = 'r', label = r'$\mathcal{B}_\kappa(\ell_h, \ell_h, \ell_s = 100)$')
plot(ll_hard, diagonal(cov_ssc_tng)/diagonal(cov_ssc_dmo)            , linewidth = 1.5, linestyle = 'solid', c = 'm', label = r'${\rm Cov}^{SSC}_\kappa(\ell_h, \ell_h)$')
# Cosmetics
fill_between(ll_hard, 0.99*ones(len(ll_hard)), 1.01*ones(len(ll_hard)), color = 'grey', alpha = 0.50)
axvline( ll_hard[where(abs(C_l_tng/C_l_dmo                             - 0.99)==min(abs(C_l_tng/C_l_dmo                             - 0.99)))], ymax = 0.48, linestyle = 'dashed', linewidth = 1., c = 'b' )
axvline( ll_hard[where(abs(B_lll_tng[:,0]/B_lll_dmo[:,0]               - 0.99)==min(abs(B_lll_tng[:,0]/B_lll_dmo[:,0]               - 0.99)))], ymax = 0.48, linestyle = 'dashed', linewidth = 1., c = 'g' )
axvline( ll_hard[where(abs(B_lll_tng[:,2]/B_lll_dmo[:,2]               - 0.99)==min(abs(B_lll_tng[:,2]/B_lll_dmo[:,2]               - 0.99)))], ymax = 0.48, linestyle = 'dashed', linewidth = 1., c = 'r' )
axvline( ll_hard[where(abs(diagonal(cov_ssc_tng)/diagonal(cov_ssc_dmo) - 0.99)==min(abs(diagonal(cov_ssc_tng)/diagonal(cov_ssc_dmo) - 0.99)))], ymax = 0.48, linestyle = 'dashed', linewidth = 1., c = 'm' )
xlabel(r'$\ell_h$' , fontsize = labelsize+2)
ylabel(r'Hydro/Gravity'        , fontsize = labelsize+2)
xscale('log')
xlim(min(ll_hard), max(ll_hard))
ylim(0.80, 1.20)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font-6.}; pylab.rcParams.update(params); legend(loc = 'upper left', ncol = 2)

fig4.savefig('fig_kappa_spectra.png')

#show()


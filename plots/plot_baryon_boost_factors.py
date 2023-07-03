import sys; sys.path.append('../'); 
from prepare_for_lenscov import *

## ========================================================
## Prepare boost factor to plot 
## ========================================================
kkhere = 10.**linspace(-2, log10(15.67724557571816), 500)

bb_z00 = zeros(len(kkhere))
bb_z05 = zeros(len(kkhere))
bb_z10 = zeros(len(kkhere))
bb_z20 = zeros(len(kkhere))
for i in range(len(kkhere)):
    bb_z00[i] = Pnl_tng_int(0.0, kkhere[i])  /  Pnl_int(0.0, kkhere[i])
    bb_z05[i] = Pnl_tng_int(0.5, kkhere[i])  /  Pnl_int(0.5, kkhere[i])
    bb_z10[i] = Pnl_tng_int(1.0, kkhere[i])  /  Pnl_int(1.0, kkhere[i])
    bb_z20[i] = Pnl_tng_int(2.0, kkhere[i])  /  Pnl_int(2.0, kkhere[i])

## ========================================================
## Plot
## ========================================================
labelsize = 26
ticksize  = 26
textsize  = 20
titlesize = 24
legend_font = 22
text_font = 24

fig1 = plt.figure(1, figsize=(8., 10.))
fig1.subplots_adjust(left=0.16, right=0.983, top=0.96, bottom=0.10, wspace = 0.28, hspace = 0.29)

# Plot hydro/gravity ratio for the power spectrum
panel = fig1.add_subplot(2,1,1)
# Add the data
plot(kkhere, bb_z00, c = 'g'         , linewidth = 1.5, linestyle = 'solid', label = r'$z=0$')
plot(kkhere, bb_z05, c = 'c'         , linewidth = 1.5, linestyle = 'solid', label = r'$z=0.5$')
plot(kkhere, bb_z10, c = 'darkorange', linewidth = 1.5, linestyle = 'solid', label = r'$z=1$')
plot(kkhere, bb_z20, c = 'k'         , linewidth = 1.5, linestyle = 'solid', label = r'$z=2$')
# Cosmetcis
xlabel(r'$k\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'$\mathcal{I}^{{\rm baryon}}(k,z)$'        , fontsize = labelsize)
xscale('log')
xlim(min(kkhere), max(kkhere))
ylim(0.75, 1.05)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'lower left', ncol = 2)

# Plot log-derivative of the power spectrum
panel = fig1.add_subplot(2,1,2)
# Add the data
kplot = 10.**linspace(-2, log10(15.67724557571816), 500)
z_now = 0.5
toplot_dmo = zeros(len(kplot))
toplot_hyd = zeros(len(kplot))
for i in range(len(kplot)):
    toplot_dmo[i] = kplot[i]*dPnl_int(z_now, kplot[i])     / Pnl_int(z_now, kplot[i])
    toplot_hyd[i] = kplot[i]*dPnl_tng_int(z_now, kplot[i]) / Pnl_tng_int(z_now, kplot[i])
plot(kplot, toplot_dmo, c = 'b', linewidth = 1.5, linestyle = 'solid', label = r'Gravity')
plot(kplot, toplot_hyd, c = 'r', linewidth = 1.5, linestyle = 'solid', label = r'Hydro')
# Cosmetics
annotate(r"$z = 0.5$"                  , xy = (0.10,0.10), xycoords='axes fraction', color = 'k', fontsize = text_font)
xlabel(r'$k\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
ylabel(r'${\rm dln} P(k) / {\rm dln k}$' , fontsize = labelsize)
xscale('log')
xlim(min(kplot), max(kplot))
ylim(-2.5, 0.5)
xticks(size = ticksize)
yticks(size = ticksize)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1)

fig1.savefig('fig_baryon_boost_factors.png')

#show()

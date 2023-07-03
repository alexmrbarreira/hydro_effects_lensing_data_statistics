import sys; sys.path.append('../');
from prepare_for_lenscov import *

docoarse = True
#docoarse = False

# ============================================================= 
# Load responses to plot
# =============================================================

# k array (same for both, so pick any)
kk = loadtxt('../lookup_tables/Resp_G1_TNG_hydro_total.dat', skiprows = 1)[:,0]

# Hydro
file_G1_hydro = loadtxt('../lookup_tables/Resp_G1_TNG_hydro_total.dat', skiprows = 1)
G1_tng300_2_hydro_total_z00 = file_G1_hydro[:,1]
G1_tng300_2_hydro_total_z05 = file_G1_hydro[:,2]
G1_tng300_2_hydro_total_z10 = file_G1_hydro[:,3]
G1_tng300_2_hydro_total_z20 = file_G1_hydro[:,4]
G1_tng300_2_hydro_total_z30 = file_G1_hydro[:,5]
G1_tng300_2_hydro_total     = [G1_tng300_2_hydro_total_z00, G1_tng300_2_hydro_total_z05, G1_tng300_2_hydro_total_z10, G1_tng300_2_hydro_total_z20, G1_tng300_2_hydro_total_z30]

# Gravity
file_G1_dmo = loadtxt('../lookup_tables/Resp_G1_TNG_dmo_total.dat', skiprows = 1)
G1_tng300_2_dmo_total_z00 = file_G1_dmo[:,1]
G1_tng300_2_dmo_total_z05 = file_G1_dmo[:,2]
G1_tng300_2_dmo_total_z10 = file_G1_dmo[:,3]
G1_tng300_2_dmo_total_z20 = file_G1_dmo[:,4]
G1_tng300_2_dmo_total_z30 = file_G1_dmo[:,5]
G1_tng300_2_dmo_total     = [G1_tng300_2_dmo_total_z00, G1_tng300_2_dmo_total_z05, G1_tng300_2_dmo_total_z10, G1_tng300_2_dmo_total_z20, G1_tng300_2_dmo_total_z30]

Nz = len(G1_tng300_2_hydro_total)

# ============================================================= 
# Interpolate to coarser array for smoother plots
# =============================================================

if(docoarse):
    kk_coarse = 10.**linspace(log10(min(kk)+1.0e-6), log10(max(kk)-1.0e-6), 30)

    # Hydro
    for i in range(Nz):
        G1_tng300_2_hydro_total[i] = interpolate.interp1d(kk, G1_tng300_2_hydro_total[i])(kk_coarse)

    # Gravity
    for i in range(Nz):
        G1_tng300_2_dmo_total[i] = interpolate.interp1d(kk, G1_tng300_2_dmo_total[i])(kk_coarse)

else:
    kk_coarse = kk

# ============================================================= 
# Load Wagner et al result
# =============================================================

file_wagner = loadtxt('../lookup_tables/Resp_G1_WagnerSims.dat', skiprows = 1)

kk_wagner     = file_wagner[:,0]
G1_wagner_z00 = file_wagner[:,1]
G1_wagner_z05 = file_wagner[:,2]
G1_wagner_z10 = file_wagner[:,3]
G1_wagner_z20 = file_wagner[:,4]
G1_wagner_z30 = file_wagner[:,5]
G1_wagner     = [G1_wagner_z00, G1_wagner_z05, G1_wagner_z10, G1_wagner_z20, G1_wagner_z30]

# ============================================================= 
# Make plot
# =============================================================

labelsize   = 26
ticksize    = 26
textsize    = 20
titlesize   = 20
text_font   = 24
legend_font = 22
alpha_c     = 0.2

x_min = min(kk_coarse)
x_max = max(kk_coarse)
y_min = -0.5
y_max = 2.5

z_names = [r"$z = 0$", r"$z = 0.5$", r"$z = 1$", r"$z = 2$", r"$z = 3$"]

fig0 = plt.figure(0, figsize=(17., 10.))
fig0.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.10, wspace = 0.35, hspace = 0.25)

for i in range(Nz):
    panel = fig0.add_subplot(2,3,1+i)
    # Add to plot
    plot(kk_coarse, G1_tng300_2_dmo_total[i]  , c = 'b', linewidth = 1.5, linestyle = 'solid', label = r'Gravity')
    plot(kk_coarse, G1_tng300_2_hydro_total[i], c = 'r', linewidth = 1.5, linestyle = 'solid', label = r'Hydro')
    plot(kk_wagner, G1_wagner[i], c = 'k', linewidth = 1.5, linestyle = 'dashed', label = r'Wagner et al (2015)')
    axhline(26./21, c = 'k', linewidth = 1.5, linestyle = 'dotted')
    # Cosmetics
    annotate(z_names[i], xy = (0.10,0.10), xycoords='axes fraction', color = 'k', fontsize = text_font)
    xlabel(r'$k\ \left[h/{\rm Mpc}\right]$'        , fontsize = labelsize)
    ylabel(r'$G_1(k)$'        , fontsize = labelsize)
    xscale('log')
    xlim(x_min, x_max)
    ylim(y_min, y_max)
    xticks(size = ticksize)
    yticks(size = ticksize)
    if(i==4):
        annotate(r"TNG300-2 resolution" , xy = (1.20,0.16), xycoords='axes fraction', color = 'k', fontsize = text_font+10)
        params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'upper right', ncol = 1, bbox_to_anchor = (2.2,0.90))

fig0.savefig('fig_G1_tng300_2.png')

#show()

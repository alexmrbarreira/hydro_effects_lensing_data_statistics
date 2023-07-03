from numpy import *
from scipy import *
from pylab import *
from scipy import interpolate, integrate, special
from scipy.interpolate import RegularGridInterpolator, UnivariateSpline
from numpy import random
import pylab,  time
import os
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
rcParams.update({'text.usetex': True, 'mathtext.fontset': 'stix'})

################################################rcParams.update({'text.usetex': False, 'mathtext.fontset': 'stix'})

def jell(ell, x): #spherical Bessel function 
    return special.spherical_jn(int(ell), x, derivative=False)

lpath = '/home/barreira/a.Other/d-github_workrepos/c-hydro_effects_lensing_data_statistics/' 

# ==================================================== #
# Compute z(chi) relation to use in interpolations
# ==================================================== #
print ('Computing z(chi) relation ... ')
Om0 = 0.3089
h   = 0.6774
c_light = 299792.458 #km/s

if (Om0 != 0.3089):
    print ('')
    print ('Om0 =', Om0, 'but should be equal to 0.3089 to match cosmology of tables in lookup_tables/')
    print ('Should you wish to modify the cosmology, you need to regenerate the spectra tables in lookup_tables/')
    print ('')
    quit()
elif (h!= 0.6774):
    print ('')
    print ('h =', h, 'but should be equal to 0.6774 to match cosmology of tables in lookup_tables/')
    print ('Should you wish to modify the cosmology, you need to regenerate the spectra tables in lookup_tables/')
    print ('')
    quit()

radian_to_arcmin = 180.*60./pi          #Multiply by this amount to convert radians to arcmin
radian_to_arcsec = radian_to_arcmin*60. #Multiply by this amount to convert radians to arcsec

z_max      = 3.0 # some sufficiently high value such that we never get to these source redshifts
zz_array = linspace(0., z_max, 10000)

def E_lcdm(a, Om0): # This is E(a) = H(a)/H0
    return sqrt(Om0*a**(-3.) + (1-Om0))
def d_ang(z_i, z_f, Nz, E, Om0): # This gives the angular diameter distance in Mpc/h between z_i and z_f for given E = H(a)/H0
    array_int_z = linspace(z_i, z_f, Nz)
    d_ang_int   = 1./E(1./(1 + array_int_z), Om0)
    return (c_light/100)*integrate.trapz(d_ang_int, array_int_z)/(1 + z_f)

chiofz = zeros(len(zz_array))
for i in range(len(zz_array)):
    chiofz[i] = d_ang(0.0, zz_array[i], 1000, E_lcdm, Om0)*(1. + zz_array[i])

chiofz_int = interpolate.interp1d(sort(zz_array), sort(chiofz))
zofchi_int = interpolate.interp1d(sort(chiofz), sort(zz_array))

# ============================================================== #
# Specify lensing kernel for desired redshift
# ============================================================== #
z_s           = 1.0
chi_s         = chiofz_int(z_s)
chi_array     = linspace(0.01, chi_s, 1000)
z_array       = zofchi_int(chi_array)
a_array       = 1./(1.+z_array)
gkernel_array = (3.*(100.)**2*Om0/2./c_light**2.) * (chi_s - chi_array)*chi_array / chi_s / a_array
gkernel_int   = interpolate.interp1d(chi_array, gkernel_array)

# ============================================================== #
# Specify ell arrays
# ============================================================== #
lbinedges    = 10.**linspace(log10(20.), log10(5000.), 51)
ll_hard      = (lbinedges[:-1] + lbinedges[1:])/2.
deltal_array = lbinedges[1:] - lbinedges[:-1]

ll_soft = array([ll_hard[0], ll_hard[11], ll_hard[14]])

print ('ll_hard = ', ll_hard)
print ('ll_soft = ', ll_soft)

# ============================================================== #
# Covariance parameters 
# ============================================================== #
Omega_S       = 0.363610*4.*pi

sigma_eps = 0.37
sounumden = 30.*radian_to_arcmin**2. #in /radian^2

# ================================================================================ #
# Load spectra and growth-only response tables and create 2D interpolators
# ================================================================================ #
print ('Making interpolators ... ')
# Create the interpolators that are needed; it does trick of setting to zero (or linear result for responses) for k = 0 (to prevent crashes/problems when k<k_min)
# This needs to be adapted if the format of the files changes
def make_interpolator_forP(path1, path2, value_at_k0):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    ktmp     = append(0., filetmp[:,0])
    Ptmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    return interpolate.interp2d(ztmp, ktmp, Ptmp, kind='linear')

def make_interpolator_forG(path1, path2, value_at_k0, index_highk, asymptote, doextrapolation):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    # Set to linear theory for low-k
    ktmp     = append(0., filetmp[:,0])
    Gtmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    if(doextrapolation):
        # Set correction for high-k
        khig = 10.**linspace(log10(max(ktmp+0.000001)), 3, 500)
        kout = append(ktmp, khig)
        Ghig = zeros([len(kout), len(ztmp)])
        Ghig[0:len(ktmp), :] = Gtmp
        for i in range(len(ztmp)):
            G_at_kmax = Ghig[len(ktmp)-1, i]
            kmax      = ktmp[-1]
            Ghig[len(ktmp)::, i] = (kmax/khig)**index_highk * (G_at_kmax - asymptote) + asymptote
        return interpolate.interp2d(ztmp, kout, Ghig, kind='cubic')
    else:
        return interpolate.interp2d(ztmp, ktmp, Gtmp, kind='cubic')

# P_L(z, k) w/ Camb
Plin_int   = make_interpolator_forP(lpath+'lookup_tables/P_L_camb.dat', lpath+'lookup_tables/P_L_camb_zvalues.dat', 0.0)

# Select if to use Halofit (0) or Emulator (1)
which_nonlinear = 1

if(which_nonlinear == 0): #Use Halofit
    print ('Using Halofit ... ')
    # P_m(z, k) , dP_m/dk, d^2P_L/dk^2, w/ Camb Halofit
    Pnl_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_camb.dat'  , lpath+'lookup_tables/P_m_camb_zvalues.dat', 0.0)
    dPnl_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_camb.dat' , lpath+'lookup_tables/P_m_camb_zvalues.dat', 0.0)
    ddPnl_int = make_interpolator_forP(lpath+'lookup_tables/ddP_m_camb.dat', lpath+'lookup_tables/P_m_camb_zvalues.dat', 0.0)
    # P_m(z, k) , dP_m/dk, d^2P_L/dk^2, w/ Camb Halofit and TNG300-2 baryon boost factor
    Pnl_tng_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_camb_TNG300_2_hydro_total.dat'  , lpath+'lookup_tables/P_m_camb_zvalues.dat', 0.0)
    dPnl_tng_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_camb_TNG300_2_hydro_total.dat' , lpath+'lookup_tables/P_m_camb_zvalues.dat', 0.0)
    ddPnl_tng_int = make_interpolator_forP(lpath+'lookup_tables/ddP_m_camb_TNG300_2_hydro_total.dat', lpath+'lookup_tables/P_m_camb_zvalues.dat', 0.0)
elif(which_nonlinear == 1): #Use Emulator
    print ('Using CosmicEmu ... ')
    # P_m(z, k) , dP_m/dk, d^2P_L/dk^2, w/ CosmicEmu
    Pnl_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_emu.dat'  , lpath+'lookup_tables/P_m_emu_zvalues.dat', 0.0)
    dPnl_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_emu.dat' , lpath+'lookup_tables/P_m_emu_zvalues.dat', 0.0)
    ddPnl_int = make_interpolator_forP(lpath+'lookup_tables/ddP_m_emu.dat', lpath+'lookup_tables/P_m_emu_zvalues.dat', 0.0)
    # P_m(z, k) , dP_m/dk, d^2P_L/dk^2, w/ CosmicEmu and TNG300-2 baryon boost factor
    Pnl_tng_int   = make_interpolator_forP(lpath+'lookup_tables/P_m_emu_TNG300_2_hydro_total.dat'  , lpath+'lookup_tables/P_m_emu_zvalues.dat', 0.0)
    dPnl_tng_int  = make_interpolator_forP(lpath+'lookup_tables/dP_m_emu_TNG300_2_hydro_total.dat' , lpath+'lookup_tables/P_m_emu_zvalues.dat', 0.0)
    ddPnl_tng_int = make_interpolator_forP(lpath+'lookup_tables/ddP_m_emu_TNG300_2_hydro_total.dat', lpath+'lookup_tables/P_m_emu_zvalues.dat', 0.0)
else:
    print (' ')
    print ('Hey, which_nonlinear not allowed !!')
    print (' ')

# G_1(z, k) and G_K(z, K) from Wagner et al and Schmidt et al simulations
G1_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_G1_WagnerSims.dat' , lpath+'lookup_tables/Resp_zvalues_G1_WagnerSims.dat' , 26./21, 0.5, -3./4, True)
GK_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_GK_SchmidtSims.dat', lpath+'lookup_tables/Resp_zvalues_GK_SchmidtSims.dat',   8./7, 0.5, -9./4, True)

# G_1(z, k) from TNG300-2; gravity-only and hydro runs 
G1_tng_dmo_int   = make_interpolator_forG(lpath+'lookup_tables/Resp_G1_TNG_dmo_total.dat'   , lpath+'lookup_tables/Resp_zvalues_G1_TNG_dmo_total.dat'   , 26./21, 0.5, -3./4, True)
G1_tng_hydro_int = make_interpolator_forG(lpath+'lookup_tables/Resp_G1_TNG_hydro_total.dat' , lpath+'lookup_tables/Resp_zvalues_G1_TNG_hydro_total.dat' , 26./21, 0.5, -3./4, True)

# ================================================================================ #
# Define response coefficients
# ================================================================================ #

# R_1(z, k) and R_K(z, k) w/ gravity nonlinear spectra, Wagner G1 and Schmidt GK.
def R_1_int(z, k):
    return 1. - (1./3)*k*dPnl_int(z, k)/Pnl_int(z, k) + G1_int(z, k)
def R_K_int(z, k):
    return GK_int(z, k) - k*dPnl_int(z, k)/Pnl_int(z, k) 

# R_1(z, k) w/ gravity nonlinear spectra, TNG dmo G1
def R_1_tng_dmo_int(z, k):
    return 1. - (1./3)*k*dPnl_int(z, k)/Pnl_int(z, k) + G1_tng_dmo_int(z, k)

# R_1(z, k) and R_K(z, k) w/ gravity nonlinear spectra with TNG baryon boost , TNG hydro G1 and Schmidt GK
def R_1_tng_int(z, k):
    return 1. - (1./3)*k*dPnl_tng_int(z, k)/Pnl_tng_int(z, k) + G1_tng_hydro_int(z, k)
def R_K_tng_int(z, k):
    return GK_int(z, k) - k*dPnl_tng_int(z, k)/Pnl_tng_int(z, k)

# ================================================================================ #
# Define function that symmetrizes matrices 
# ================================================================================ #
# Function to symmetrize calculation of the covariance matrices. For instanc can speed up things by computing only the triangle k1>=k2
# Useufl for squeezed covraince because the matrix has to be symmetric wrt k1=k2, but things depends on whether k1><k2; we compute the matrix assuming k1>k2, but then symmetrize it 
def symmetrize_matrix(matrix):
    if(len(matrix[:,0]) != len(matrix[0,:])):
        print ('Not square matrix!'); quit()
    n = len(matrix[:,0])
    new_matrix = copy(matrix)
    for j in range(n):
        for i in range(1+j, n):
            new_matrix[i, j] = matrix[j, i]
    return new_matrix


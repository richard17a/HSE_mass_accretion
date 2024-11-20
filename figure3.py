"""
     author: Maylis Landeau (IPGP)
"""

import numpy as np    
import matplotlib.pyplot as plt
import matplotlib
import tol_colors as tc
cset = tc.tol_cset('muted')
import os,sys
plt.close("all")
from utils import set_size


matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#Colorblind-friendly set
#https://personal.sron.nl/~pault/
Color_1 = cset[2]
Color_2 = cset[7]
Color_3 = cset[3]
Color_4 = cset[1]
Alpha_1 = 0.2

# -------------------------------
# -------------------------------
# Input parameters
# -------------------------------
# -------------------------------

Fontsize = 10 
Fontsize_small = 8.6
k =   1.e-6 #Thermal diffusivity, m^2/s
Cp = 1.e3 # Heat capacity, in J / kg / K
H_all = [2e6,1.e4] # Magma ocean thickness, in m
H = 2.e6
VMO = 1./2. * 4./3.*np.pi * H**3 # Magma ocean volume in m^3, hemispheric magma pond

Rp = 6.3e6 # Planet radius in m 
Rc = 3.48e6 #Core radius in m
#VMO = 4./3.*np.pi*(Rp**3-(Rp-H)**3) # Magma ocean volume in m^3
Vmantle = 4./3.*np.pi*(Rp**3-Rc**3) # Mantle volume in m^3
g = 9.8 #gravity, m/s^2
alpha = 3e-5 #mantle thermal expansion coefficient, in 1/K
mu_fixed = 0.05 #Dynamic viscosity, Pa s
rho_s = 4500.0 #Mantle density, kg/m^3
nu_fixed = mu_fixed / rho_s  #Mantle kinematic viscosity, m^2/s
rho_m = 9000.0 # metal density, kg/m^3
K = k * rho_s * Cp # thermal condcutivity
F_fixed = [0.3e6,1.e6,3.e6] #Heat flux, W / m^2
Color_F = [Color_1,Color_2,Color_3]
Mass_metal_needed_Earth = 1.e22 #Mass of metal needed to explain Earth HSE (kg)
Mass_metal_needed_Earth_min = 0.6e22 #Mass of metal needed to explain Earth HSE (kg)
Mass_metal_needed_Earth_max = 1.3e22 #Mass of metal needed to explain Earth HSE (kg)
Phi_Earth_min = Mass_metal_needed_Earth_min / (rho_m * Vmantle)
Phi_Earth_max = Mass_metal_needed_Earth_max / (rho_m * Vmantle)
Cn = 0.3

# Coefficient in volume fraction in suspension, Solomatov et al. 1993
eps_all = [0.2e-2,0.9e-2] 
eps_Lavorel = 0.25
Rac = 1.e3 #Critical Rayleigh number


# -------------------------------
# -------------------------------
# Scalings and functions
# -------------------------------
# -------------------------------


def Entrained_Fraction_Drag(eps,Up,alpha,F,Cp,rho_m,rho_s) : 
    #Solomatov & Stevenson 1993, with terminal velocity
    Phi = eps*alpha*F/(Cp*(rho_m-rho_s)*Up)
    
    return Phi



def Terminal_Velocity(rho_m,rho_s,g,r,Cd): 
    #r = drop radius
    Up = (8./3.*(rho_m-rho_s)/rho_s*g*r/Cd)**0.5
    return Up


def Drag_coeff (r,mu,Cn,g,rho_m,rho_s) :
    # Samuel 2012, modifying Cn
    Ustokes = 2./9.*(g*(rho_m-rho_s)*r**2/mu)
    nu = mu / rho_s
    Rep = Ustokes * r/nu
    Cd = 12./Rep + Cn
    return Cd

def Stokes_velocity (r,mu,g,rho_m,rho_s) :
    # Samuel 2012, modifying Cn
    Ustokes = 2./9.*(g*(rho_m-rho_s)*r**2/mu)
    return Ustokes


# -------------------------------
# -------------------------------
# Plotting
# -------------------------------
# -------------------------------

fig_width, fig_height = set_size('thesis', 1, (1, 1))

fig2, ax2 = plt.subplots(1,2,figsize=(1.5 * fig_width, fig_height))

for i_eps in range(len(eps_all)) : 
    eps = eps_all[i_eps]
    log_R =  np.arange(-6,-1,0.01) 
    R = 10**(log_R) # Drop size in m
    
    plt.subplot(1,2,i_eps+1)

    #Volume fraction in suspension, Solomatov et al. 1993, fixed flux
    for i_F in range(len(F_fixed)): 
        
        
        Cd =  Drag_coeff (R,mu_fixed,Cn,g,rho_m,rho_s) 
        Up = Terminal_Velocity(rho_m,rho_s,g,R,Cd)
        Phi_shields_F = Entrained_Fraction_Drag(eps,Up,alpha,F_fixed[i_F],Cp,rho_m,rho_s) 
        Vsuspended_shields_F = VMO * Phi_shields_F
        Mass_suspended_shields_F = rho_m * Vsuspended_shields_F
        
    
        if F_fixed[i_F] < 1.e6 : 
            plt.loglog(2.*R,Phi_shields_F,color = Color_F[i_F],label=r'$F='+ str(int(F_fixed[i_F]/1.e5)) + r'\cdot 10^5$ ${\rm W\,m}^{-2}$')
            
        else : 
            plt.loglog(2.*R,Phi_shields_F,color = Color_F[i_F],label=r'$F='+ str(int(F_fixed[i_F]/1.e6)) + r'\cdot 10^6$ ${\rm W\,m}^{-2}$')
    
    plt.axhspan(Phi_Earth_min,Phi_Earth_max,0. ,1., color = Color_3, alpha=Alpha_1, lw=0)
    plt.legend(fontsize=12, loc=3,framealpha=1,edgecolor='w')
    plt.xlabel(r'Drop diameter [m]',fontsize=16)
    if i_eps == 0:
        plt.ylabel(r'Suspended metal fraction',fontsize=16)
    plt.xlim((1.e-6,1e-1))
    plt.ylim((1.e-9,0.003))
    plt.text(3.e-3, 0.00085, r"Earth's HSEs",fontsize=12,rotation=0, color = Color_3)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=0.5)

plt.tight_layout()

ttt = ['a', 'b']
axes = fig2.get_axes()
for p,l in zip(axes, ttt):
    p.annotate(l, xy=(-0., 1.04), xycoords="axes fraction", fontsize=Fontsize, weight = 'bold')
    p.tick_params(axis='both', which='major', labelsize=13)

# plt.savefig('./figures/figure2.pdf', bbox_inches='tight', format='pdf') 

plt.show()

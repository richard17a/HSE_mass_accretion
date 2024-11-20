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


def Entrained_Fraction_Drag(eps,Up,alpha,F,Cp,rho_m,rho_s) : 
    #Solomatov & Stevenson 1993, with terminal velocity
    
    Phi = eps*alpha*F/(Cp*(rho_m-rho_s)*Up)
    
    return Phi


def Terminal_Velocity(rho_m,rho_s,g,r,Cd): 
    #d = drop radius
    
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


#Functions for turbulent convection
#----------------------------------


def HeatFlux_Convection(alpha,g,L,rho_s,Cp,k,K,nu,DT) :
    #Compute convective speed in magma ocean
    # Scaling from Solomatov 2000 (and reference within), as in Lichtenberg 2021 for comparison
    b_conv = 0.089
    Ra = alpha*rho_s*g*DT*L**3/(k*nu)
    F = b_conv * K * DT * Ra**(1./3.)/L
    
    return F


def Velocity_Convection_soft(alpha,g,L,rho_s,Cp,k,K,nu,F) :
    #Compute convective speed in magma ocean
    # Scaling from Solomatov 2000 (and reference within), as in Lichtenberg 2021 for comparison
    a_conv = 0.6
    U = a_conv * (alpha*g*L*F/(rho_s*Cp))**(1./3.)
    
    return U


def Velocity_Convection_hard(alpha,g,L,rho_s,Cp,k,K,mu,F) :
    #Compute convective speed in magma ocean
    # Scaling from Solomatov 2000, 2015 (and reference within), as in Lichtenberg 2021 for comparison
    x = 2.5*np.log(0.086*rho_s*H/mu* (alpha*g*L*F/(rho_s*Cp))**(1./3.)) + 6
    a_hard = 0.086*x
    U = a_hard * (alpha*g*L*F/(rho_s*Cp))**(1./3.)
    
    return U


def HeatFlux_Velocity_Convection(alpha,g,L,rho_s,Cp,k,K,nu,U) :
    #Compute convective speed in magma ocean
    # Scaling from Solomatov 2000 (and reference within), as in Lichtenberg 2021 for comparison
    a_conv = 0.6
    F = (U/a_conv)**3 / alpha / g / L * rho_s * Cp  
    
    return F

def x_Friction(U,H,mu_fixed,rho) : 
    #Solomatov, 2015, TOG
    x_parameter = np.arange(1.,1.e3,1.e-2)
    x_model = 2.5*np.log(rho*H*U/(mu_fixed)/x_parameter)+6
    X = np.where(x_model<x_parameter)
    x = x_parameter[X[0][0]]

    return x


def Stresses_Reynolds(U,H,mu_fixed,rho) : 
    #Solomatov, 2015, TOG
    x = x_Friction(U,H,mu_fixed,rho) 
    u = U / x
    tau = rho*u**2
    
    return tau


def Stresses_Buoyancy(F,mu_fixed,rho,alpha,g,Cp) : 
    #Solomatov et al. 1993
    tau = (mu_fixed*alpha*g*F/Cp)**0.5
    
    return tau

# -------------------------------
# -------------------------------
# Input parameters
# -------------------------------
# -------------------------------

Fontsize = 10 
Fontsize_small = 9 
k =   1.e-6 #Thermal diffusivity, m^2/s
Cp = 1.e3 # Heat capacity, in J / kg / K
H = 3e4 # Magma ocean thickness, in m
Rp = 6.3e6 # Planet radius in m 
VMO = 4./3.*np.pi*(Rp**3-(Rp-H)**3) # Magma ocean volume in m^3
DT = 500 # Temperature jump across the mantle, in K 
g = 10 #gravity, m/s^2
alpha = 3e-5 #mantle thermal expansion coefficient, in 1/K
log_mu =  np.arange(-2,0,0.01) #log10(Mdynamic viscosity)
mu = 10**(log_mu) #Dynamic viscosity, Pa s
mu_fixed = 0.02 #Dynamic viscosity, Pa s
rho_s = 4500.0 #Mantle density, kg/m^3
nu = mu/rho_s #Mantle kinematic viscosity, m^2/s
nu_fixed = mu_fixed / rho_s  #Mantle kinematic viscosity, m^2/s
rho_m = 9000.0 # metal density, kg/m^3
K = k * rho_s * Cp # thermal condcutivity
Shields_c = [0.1,0.2] # Critical Shields number, criterion for re-entrainment
Shields_c = 0.15 # Critical Shields number, criterion for re-entrainment
U_free =  np.arange(1,40,0.1) #Convective speed, m/s
U_fixed = [1,3.,10] #Velocity, m/s
log_F_free = np.arange(5,7,0.1)
F_free =10**log_F_free  #Heat flux, W / m^2
log_H_free =  np.arange(4,6.8,0.1)
H_free = 10**log_H_free
F_fixed = [3.e5,1.e6,3.e6] #Heat flux, W / m^2
Color_F = [Color_1,Color_2,Color_3]

Mass_metal_needed_Earth = 1.e22 #Mass of metal needed to explain Earth HSE (kg)
Mass_metal_needed_Earth_min = 0.6e22 #Mass of metal needed to explain Earth HSE (kg)
Mass_metal_needed_Earth_max = 1.3e22 #Mass of metal needed to explain Earth HSE (kg)
Cn = 0.3

# Coefficient in volume fraction in suspension, Solomatov et al. 1993
eps = 0.3e-2
eps_Lavorel = 0.25
Rac = 1.e3 #Critical Rayleigh number

# -------------------------------
# -------------------------------
# Plotting
# -------------------------------
# -------------------------------

fig_width, fig_height = set_size('thesis', 1, (1, 1))

fig, ax = plt.subplots(1,1,figsize=(fig_width, fig_height))
tau = H_free*0
for i_F in range(len(F_fixed)) : 
    F = F_fixed[i_F]
    for i_H in range(len(H_free)) : 
        H = H_free[i_H]
        U_soft =  Velocity_Convection_soft(alpha,g,H,rho_s,Cp,k,K,nu_fixed,F)
        U_hard = Velocity_Convection_hard(alpha,g,H,rho_s,Cp,k,K,mu_fixed,F)
        tau_R_i_soft = Stresses_Reynolds(U_soft,H,mu_fixed,rho_s) 
        tau_R_i_hard = Stresses_Reynolds(U_hard,H,mu_fixed,rho_s) 
        tau_B_i = Stresses_Buoyancy(F,mu_fixed,rho_s,alpha,g,Cp) 
        tau_i = np.max((tau_R_i_hard,tau_B_i))
        tau[i_H] = tau_i
        D_c = tau/(rho_m-rho_s)/g/Shields_c
        
    if F < 1.e6 : 
        plt.loglog(D_c,H_free,'-',color=Color_F[i_F],label=r'$F='+ str(int(F_fixed[i_F]/1.e5)) + r'\cdot 10^5~{\rm W\,m}^{-2}$')

        
    else : 
        plt.loglog(D_c,H_free,'-',color=Color_F[i_F],label=r'$F='+ str(int(F_fixed[i_F]/1.e6)) + r'\cdot 10^6~{\rm W\,m}^{-2}$')

    plt.ylabel(r'Magma pond depth [m]', fontsize=13)
    plt.xlabel('Drop diameter [m]', fontsize=13)
    
plt.legend(frameon=False, fontsize=11)
plt.ylim((1.e4,3.e6))
plt.xlim((1.e-4,1))

plt.gca().tick_params(axis='both', which='major', labelsize=12)

plt.text(2.e-2, 1.3e4, r'Settling',fontsize=12,rotation=55)
plt.text(7.e-4,1.3e4,  r'Re-entrainment',fontsize=12,rotation=55)

# plt.savefig(figure_name, bbox_inches = 'tight' )

plt.show()


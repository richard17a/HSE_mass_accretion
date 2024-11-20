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
from cmcrameri import cm


matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# -------------------------------
# -------------------------------
# Input parameters
# -------------------------------
# -------------------------------

c = [0.05,0.2] #Prefactor in velocity scaling, from 0.05 to 0.2
k = 1.e-6 #Thermal diffusivity, m^2/s
H = 2.9e6 # Solid mantle thickness, in m
DT = [500,2500] # Temperature jump across the mantle, in K 
g = 10 #gravity, m/s^2
alpha = 3e-5 #mantle thermal expansion coefficient, in 1/K

log_mu =  np.arange(15,20,0.01) #log10(Mdynamic viscosity)
log_mu_mantle = np.array([16.0,19.0])#log10(Mantle dynamic viscosity)
mu = 10**(log_mu) #Dynamic viscosity, Pa s
mu_mantle = 10**(log_mu_mantle) #Mantle dynamic viscosity, Pa s
rho = 4500.0 #Mantle density, kg/m^3
nu = mu/rho #Mantle kinematic viscosity, m^2/s
rho_m = 9000.0 # metal density, kg/m^3
Shields_c = [0.2,0.4] # Critical Shields number, criterion for re-entrainment
# Coefficient in volume fraction in suspension, Solomatov et al. 1993
eps = 0.3e-2
Mass_metal_needed_Earth = 1.e22 #Mass of metal needed to explain Earth HSE (kg)

Rp = 6.3e6 # Planet radius in m 
VMantle = 4./3.*np.pi*(Rp**3-(Rp-H)**3) # Magma ocean volume in m^3


# -------------------------------
# -------------------------------
# Scalings
# -------------------------------
# -------------------------------

Ra1 = (alpha*g*DT[0]*H**3)/(nu*k)
Ra2 = (alpha*g*DT[1]*H**3)/(nu*k)

U1 =  c[0]*k/H* Ra1**(2./3.)
U2 =  c[1]*k/H* Ra2**(2./3.)

R_c1 = mu * U1 / ((rho_m-rho)* g * H * Shields_c[1]) #Critical diapir radius for re-entrainment
R_c2 = mu * U2 / ((rho_m-rho)* g * H * Shields_c[0]) #Critical diapir radius for re-entrainment

#Volume fraction in suspension, Solomatov et al. 1993
R = 10 # Drop size in m
tau1 = mu * U1 / H
tau2 = mu * U2 / H
Phi_shields1 = 18 * eps * (tau1/((rho_m-rho)* g * R))**2
Phi_shields2 = 18 * eps * (tau2/((rho_m-rho)* g * R))**2

#Equivalent spherical radius
Vsuspended_shields1 = VMantle  * Phi_shields1
Vsuspended_shields2 = VMantle  * Phi_shields2
Mass_suspended_shields1 = rho_m * Vsuspended_shields1
Mass_suspended_shields2 = rho_m * Vsuspended_shields2


# -------------------------------
# -------------------------------
# Plotting
# -------------------------------
# -------------------------------

fig_width, fig_height = set_size('thesis', 1, (1, 1))

fig, ax = plt.subplots(1,1,figsize=(fig_width, fig_height))

plt.loglog(2*R_c1,mu,'--',color=cm.bamako(0.33))
plt.loglog(2*R_c2,mu,':',color=cm.bamako(0.33))
plt.ylabel('Mantle viscosity [Pa s]', fontsize=13)
plt.xlabel('Metal diapir diameter [m]', fontsize=13)
plt.ylim((1.e15,1.e20))
plt.xlim((0.4,500))
plt.xticks([1, 10, 100], labels=[1, 10, 100])
X = np.where(mu>mu_mantle[0])
i_1 = X[0][0]
X = np.where(mu<mu_mantle[1])
i_2 = X[0][-1]

ax.fill_betweenx(mu_mantle, [2*R_c1[i_1],2*R_c1[i_2]], [2*R_c2[i_1],2*R_c2[i_2]], alpha=0.2, facecolor=cm.bamako(0.33))

plt.text(65, 5.e16, r'Settling',fontsize=12,rotation=48)
plt.text(0.8,6.e16,  r'Possible entrainment',fontsize=12,rotation=48)

# plt.savefig('./figures/figure4.pdf', bbox_inches='tight', format='pdf') 

plt.show()

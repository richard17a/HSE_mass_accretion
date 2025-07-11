import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from cmcrameri import cm

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
        width_pt = width_pt / 72.27
    elif width == 'beamer':
        width_pt = 307.28987
        width_pt = width_pt / 72.27
    else:
        width_pt = width

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2


    fig_width_in = width_pt * fraction
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)
    

fig_width, fig_height = set_size('thesis', 1, (1, 1))


def main():

    M_earth = 5.97e24
    M_moon = 7.35e22
    rho_imp = 3000

    M_HSE = 0.005 * M_earth 

    p_SI_low = 1.3
    p_SI_high = 1.6

    alpha_SI_low = 3 * p_SI_low - 2
    alpha_SI_high = 3 * p_SI_high - 2

    alpha_MAB_low = 2.1 # Tsirvoulis et al 2018
    alpha_MAB_high = 3.3 # SDSS inner MAB Pena et al 2020, Maeda et al 2021

    D_maxs = np.logspace(np.log10(5), np.log10(3e6), 1000)
    D_min = 0.1

    K = (3 * M_HSE) / (np.pi * rho_imp * (1000 ** 0.5 - D_min ** 0.5))
    M_tot_35_1000 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.025
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_0025 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.05
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_005 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.1
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_01 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.2
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_02 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.3
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_03 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.4
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_04 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.5
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_05 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.6
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_06 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.7
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_07 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.8
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_08 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 0.9
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_09 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    f0 = 1.0
    K = (3 * M_HSE / (np.pi * rho_imp)) / (f0 * (D_maxs ** 0.5) + (1 - f0) * (1000 ** 0.5) - (D_min ** 0.5))
    M_tot_10 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    K = (3 - alpha_SI_low + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_SI_low + 1) - D_min ** (3 - alpha_SI_low + 1)))
    M_tot_SI_low = K * np.pi * rho_imp * (D_maxs ** (3 - alpha_SI_low + 1) - D_min ** (3 - alpha_SI_low + 1)) / 6 / (3 - alpha_SI_low + 1)

    K = (3 - alpha_SI_high + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_SI_high + 1) - D_min ** (3 - alpha_SI_high + 1)))
    M_tot_SI_high = K * np.pi * rho_imp * (D_maxs ** (3 - alpha_SI_high + 1) - D_min ** (3 - alpha_SI_high + 1)) / 6 / (3 - alpha_SI_high + 1)

    D_max_MAB = np.logspace(np.log10(5), 5, 1000)

    K = (3 - alpha_MAB_low + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_MAB_low + 1) - D_min ** (3 - alpha_MAB_low + 1)))
    M_tot_MAB_low = K * np.pi * rho_imp * (D_max_MAB ** (3 - alpha_MAB_low + 1) - D_min ** (3 - alpha_MAB_low + 1)) / 6 / (3 - alpha_MAB_low + 1)

    K = (3 - alpha_MAB_high + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_MAB_high + 1) - D_min ** (3 - alpha_MAB_high + 1)))
    M_tot_MAB_high = K * np.pi * rho_imp * (D_max_MAB ** (3 - alpha_MAB_high + 1) - D_min ** (3 - alpha_MAB_high + 1)) / 6 / (3 - alpha_MAB_high + 1)

    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cm.bamako, norm=norm)
    sm.set_array([])

    _ = plt.figure(figsize=(fig_width, fig_height))

    f_values = [0, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    M_tot_values = [M_tot_35_1000, M_tot_0025, M_tot_005, M_tot_01, M_tot_02, M_tot_03, M_tot_04,
                    M_tot_05, M_tot_06, M_tot_07, M_tot_08, M_tot_09, M_tot_10]

    for f0, M_tot in zip(f_values, M_tot_values):
        if f0 == 0.6:
            plt.plot(D_maxs / 1e3, M_tot / M_HSE, color=cm.bamako(f0), zorder=10, lw=2, label=rf'$f_0 \approx 0.6$ (Fischer-Gödde et al., 2020)')
        elif f0 == 0.3:
            plt.plot(D_maxs / 1e3, M_tot / M_HSE, color=cm.bamako(f0), zorder=10, lw=2, label=rf'$f_0 \approx 0.6$ (Worsham & Kleine, 2021)')
        else:
            plt.plot(D_maxs / 1e3, M_tot / M_HSE, color=cm.bamako(f0), zorder=10, lw=0.75)

    plt.axvline(946, c='tab:gray', ls=':', zorder=0)
    plt.text(950, 2.7e24 / M_HSE, r"Ceres", rotation=270, fontsize=10, color='tab:gray')

    plt.axhline(M_earth / M_HSE, ls='--', c='tab:red', zorder=0)
    plt.text(0.65, 4e24 / M_HSE, r"1 Earth mass", fontsize=11, color='tab:red')

    plt.axhline(M_moon / M_HSE, ls='--', c='tab:red', zorder=0)
    plt.text(0.65, 9e22 / M_HSE, r"1 Moon mass", fontsize=11, color='tab:red')

    plt.axhspan(0.004 * M_earth / M_HSE, 0.006 * M_earth / M_HSE, color='darkgreen', alpha=0.25, zorder=1)
    plt.text(200, 0.875, r"$M_{\rm HSE, \oplus}$", fontsize=11, color='darkgreen', zorder=20)

    plt.xlim(4e-1, 3e3)
    plt.ylim(2e22 / M_HSE, 8e24 / M_HSE)

    plt.xscale('log')
    plt.yscale('log')

    plt.xticks([1, 10, 100, 1000], labels=[1, 10, 100, 1000])
    plt.yticks([1, 10, 100], labels=[1, 10, 100])

    plt.xlabel(r'$D_{\rm max}$ [km]', fontsize=13)
    plt.ylabel(r'Total accreted mass $[M_{\rm HSE, \,\oplus}]$', fontsize=13)

    cbar = plt.colorbar(sm, ax=plt.gca(), pad=0.02)
    
    cbar.set_label(r'$f_o$', fontsize=13)
    cbar.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    cbar.ax.minorticks_on()

    plt.minorticks_on()

    # plt.savefig('figures/figure7.pdf', bbox_inches='tight', format='pdf')

    plt.show()


if __name__ == "__main__":

    main()

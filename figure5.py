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

    D_crit = np.logspace(1, 4, 1000)

    K = (3 * M_HSE) / (np.pi * rho_imp * (D_crit ** 0.5 - D_min ** 0.5))
    M_tot_35_5e5 = K * np.pi * rho_imp * (5e5 ** 0.5 - D_min ** 0.5) / 3
    M_tot_35_1e6 = K * np.pi * rho_imp * (1e6 ** 0.5 - D_min ** 0.5) / 3

    K = (3 * M_HSE) / (np.pi * rho_imp * (1000 ** 0.5 - D_min ** 0.5))
    M_tot_35_1000 = K * np.pi * rho_imp * (D_maxs ** 0.5 - D_min ** 0.5) / 3

    K = (3 - alpha_SI_low + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_SI_low + 1) - D_min ** (3 - alpha_SI_low + 1)))
    M_tot_SI_low = K * np.pi * rho_imp * (D_maxs ** (3 - alpha_SI_low + 1) - D_min ** (3 - alpha_SI_low + 1)) / 6 / (3 - alpha_SI_low + 1)

    K = (3 - alpha_SI_high + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_SI_high + 1) - D_min ** (3 - alpha_SI_high + 1)))
    M_tot_SI_high = K * np.pi * rho_imp * (D_maxs ** (3 - alpha_SI_high + 1) - D_min ** (3 - alpha_SI_high + 1)) / 6 / (3 - alpha_SI_high + 1)

    D_max_MAB = np.logspace(np.log10(5), 5, 1000)

    K = (3 - alpha_MAB_low + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_MAB_low + 1) - D_min ** (3 - alpha_MAB_low + 1)))
    M_tot_MAB_low = K * np.pi * rho_imp * (D_max_MAB ** (3 - alpha_MAB_low + 1) - D_min ** (3 - alpha_MAB_low + 1)) / 6 / (3 - alpha_MAB_low + 1)

    K = (3 - alpha_MAB_high + 1) * (6 * M_HSE) / (np.pi * rho_imp * (1000 ** (3 - alpha_MAB_high + 1) - D_min ** (3 - alpha_MAB_high + 1)))
    M_tot_MAB_high = K * np.pi * rho_imp * (D_max_MAB ** (3 - alpha_MAB_high + 1) - D_min ** (3 - alpha_MAB_high + 1)) / 6 / (3 - alpha_MAB_high + 1)


    fig, ax = plt.subplots(1, 2, figsize=(1.75 * fig_width, fig_height))

    plt.subplot(1,2,1)

    plt.plot(D_maxs / 1e3, M_tot_35_1000 / M_HSE, color=cm.bamako(0.4), label=r'$D_{\rm melt} = 100\,$m', zorder=0)

    plt.fill_between(D_max_MAB / 1e3, M_tot_MAB_low / M_HSE, M_tot_MAB_high / M_HSE, color=cm.bamako(0.7), alpha=0.4, zorder=0)
    plt.fill_between(D_maxs / 1e3, M_tot_SI_low / M_HSE, M_tot_SI_high / M_HSE, color=cm.bamako(0.1), alpha=0.4, zorder=4)

    plt.text(1.3, 5, 'streaming instability', rotation=63, fontsize=13, color=cm.bamako(0.1), alpha=0.8)
    plt.text(5.0, 4.5, r'main asteroid belt', rotation=34, fontsize=13, color=cm.bamako(0.7))
    plt.text(37., 4, r'collisional cascade', rotation=25, fontsize=13, color=cm.bamako(0.4))

    plt.axvline(946, c='tab:gray', ls=':', zorder=0)
    plt.text(950, 7.5e24 / M_HSE, r"Ceres", rotation=270, fontsize=11, color='tab:gray')

    plt.axhline(M_earth / M_HSE, ls='--', c='tab:red')
    plt.text(0.55, 7e24 / M_HSE, r"1 Earth mass", fontsize=13, color='tab:red')

    plt.axhline(M_moon / M_HSE, ls='--', c='tab:red')
    plt.text(80, 9e22 / M_HSE, r"1 Moon mass", fontsize=13, color='tab:red')

    plt.axhspan(0.004 * M_earth / M_HSE, 0.006 * M_earth / M_HSE, color='darkgreen', alpha=0.25, zorder=1)
    plt.text(200, 0.875, r"$M_{\rm HSE, \oplus}$", fontsize=13, color='darkgreen')

    plt.text(4.5e-1, 650, '(a)', weight='bold')

    plt.xlim(4e-1, 3e3)
    plt.ylim(2e22 / M_HSE, 1.7e25 / M_HSE)

    plt.xscale('log')
    plt.yscale('log')

    plt.xticks([1, 10, 100, 1000], labels=[1, 10, 100, 1000])
    plt.yticks([1, 10, 100], labels=[1, 10, 100])

    plt.xlabel(r'$D_{\rm max}$ [km]', fontsize=16)
    plt.ylabel(r'Total accreted mass $[M_{\rm HSE, \,\oplus}]$', fontsize=16)

    plt.gca().tick_params(axis='both', which='major', labelsize=12)

    plt.minorticks_on()

    ylims = ax[0].get_ylim()

    ax1 = plt.gca().twinx()

    ax1.set_yscale('log')

    ax1.set_ylim([ylims[0] * M_HSE, ylims[1] * M_HSE])
    ax1.set_yticks([])

    plt.subplot(1,2,2)

    plt.plot(D_crit / 1e3, M_tot_35_5e5 / M_HSE, color=cm.bamako(0.4), ls='--', label=r'$D_{\rm max} = 500\,$km')
    plt.plot(D_crit / 1e3, M_tot_35_1e6 / M_HSE, color=cm.bamako(0.4), ls='-', label=r'$D_{\rm max} = 1000\,$km')

    plt.axhspan(0.004 * M_earth / M_HSE, 0.006 * M_earth / M_HSE, color='darkgreen', alpha=0.25, zorder=1)

    plt.text(3e0, 0.875, r"$M_{\rm HSE, \oplus}$", fontsize=13, color='darkgreen')

    plt.text(8.5e-3, 650, '(b)', weight='bold')

    plt.axhline(M_earth / M_HSE, c='tab:red', ls='--')
    plt.axhline(M_moon / M_HSE, c='tab:red', ls='--')

    plt.xlabel(r'$D_{\rm crit}$ [km]', fontsize=16)

    plt.xscale('log')
    plt.yscale('log')

    plt.minorticks_on()

    plt.legend(frameon=False, loc=(0.04, 0.2), fontsize=12)

    plt.ylim(ylims)
    plt.yticks([])

    plt.gca().tick_params(axis='both', which='major', labelsize=12)

    ax2 = plt.gca().twinx()
    ax2.set_yscale('log')

    ax2.set_ylim([ylims[0] * M_HSE, ylims[1] * M_HSE])
    ax2.set_ylabel(r'Total accreted mass [kg]', fontsize=16, rotation=270, labelpad=16)

    ax2.set_xticks([1e-2, 1e-1, 1e0, 1e1], labels=[0.01, 0.1, 1, 10])

    ax2.tick_params(axis='both', which='major', labelsize=12)

    plt.minorticks_on()

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05)

    # plt.savefig('./figures/figure5.pdf', bbox_inches='tight', format='pdf')

    plt.show()


if __name__ == "__main__":

    main()

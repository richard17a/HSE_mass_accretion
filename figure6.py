import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import ScalarFormatter
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
    """
    Docstring
    """
    AU = 1.496e11
    G = 6.67e-11
    M_e = 5.972e24
    M_m = 7.35e22
    a_m = 0.3844e9
    R_m = 1736e3
    M_sun = 1.989e30
    U_s = 4472

    v_kep = np.sqrt(G * M_e / a_m)
    v_esc_m = np.sqrt(2 * G * M_m / R_m)

    v_infty = np.logspace(-1, np.log10(2), 1000) * 11.19e3
    v_imp_m = np.sqrt(v_infty**2 + v_esc_m**2 + 3*v_kep**2)

    R = 5
    M = v_imp_m / U_s
    Fr = v_imp_m ** 2 / (1.62 * R)

    z10 = (1.092 * (1 + 0.11 * M ** 2) ** -0.25) * R * Fr ** 0.25

    R = 50
    M = v_imp_m / U_s
    Fr = v_imp_m ** 2 / (1.62 * R)

    z100 = (1.092 * (1 + 0.11 * M ** 2) ** -0.25) * R * Fr ** 0.25

    R = 500
    M = v_imp_m / U_s
    Fr = v_imp_m ** 2 / (1.62 * R)

    z1000 = (1.092 * (1 + 0.11 * M ** 2) ** -0.25) * R * Fr ** 0.25

    _ = plt.figure(figsize=(fig_width, fig_height))

    # plt.plot(v_infty / 11.19e3, z10, color=cm.bamako(0.2), label=r'$D_{\rm imp} = 10\,$m')
    plt.plot(v_infty / 11.19e3, z100, color=cm.bamako(0.5), label=r'$D_{\rm imp} = 100\,$m')
    plt.plot(v_infty / 11.19e3, z1000, color=cm.bamako(0.8), label=r'$D_{\rm imp} = 1000\,$m')

    plt.axhspan(34e3, 43e3, color='tab:gray', alpha=0.25)
    plt.text(0.6, 36e3, "Lunar crust (Wieczorek et al., 2013)", size=9)

    plt.xscale('log')
    plt.yscale('log')

    plt.ylim(5e2, )
    plt.yticks([1e3, 1e4], labels=[1, 10])
    plt.xticks([1e-1, 1e0], labels=[0.1, 1])

    plt.xlabel(r'$v_{\rm rel} / v_{\rm esc, \oplus}$', fontsize=13)
    plt.ylabel('Crater depth [km]', fontsize=13)

    plt.minorticks_on()

    plt.legend(loc=(0.01, 0.65), frameon=False)

    ax1 = plt.gca()
    ax2 = ax1.twiny()

    ax2.plot(v_imp_m / 1e3, z10, linewidth=0)

    ax2.set_xscale('log')

    ax2.set_xlabel(r'$v_{\rm imp, L}\,[{\rm km}\,{\rm s}^{-1}]$', fontsize=13)

    ax2.set_xticks([3, 4,  6, 10, 20])
    ax2.set_xticklabels(['3', '4', '6', '10', '20'])

    # plt.savefig('./figures/figure6.pdf', format='pdf', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    
    main()


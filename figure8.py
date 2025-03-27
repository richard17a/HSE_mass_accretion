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


def calc_u_i(v_imp, S_i, S_t, C_i, C_t, rho_i, rho_t):

  A = rho_i * S_i - rho_t * S_t
  B = rho_i * C_i + rho_t * C_t + 2 * rho_t * S_t * v_imp
  C = -rho_t * v_imp * (C_t + S_t * v_imp)

  u_i = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)

  return u_i


def calc_Pmax(v_imp, S_i, S_t, C_i, C_t, rho_i, rho_t):

  u_i = calc_u_i(v_imp, S_i, S_t, C_i, C_t, rho_i, rho_t)

  P_max = rho_i * u_i * (C_i + S_i * u_i)

  return P_max


def plot_vimp_Pmax():

  v_imp = np.linspace(1, 30, 1000)

  rho_i = 3314
  rho_t = 2629

  C_i = 5.43
  C_t = 3.816

  S_i = 1.34
  S_t = 1.28

  P_max = calc_Pmax(v_imp, S_i, S_t, C_i, C_t, rho_i, rho_t)

  _ = plt.figure(figsize=(fig_width, fig_height))

  plt.plot(v_imp, P_max / 1e3)

  plt.axhline(105, ls='--', color='tab:gray')
  plt.text(20, 120, 'incipient melting pressure', color='tab:gray')

  plt.yscale('log')

  plt.xlabel(r'$v_{\rm imp}$ [${\rm km}\,{\rm s}^{-1}$]', fontsize=13)
  plt.ylabel(r'$P_{\rm max}$ [GPa]', fontsize=13)

  plt.minorticks_on()

  plt.show()


def main():

    v_rel = np.linspace(0, 50, 1000)

    vimp_E = np.sqrt(11.19 ** 2 + v_rel ** 2)
    vimp_L = np.sqrt(2.38 ** 2 + 3 + 1.022 ** 2 + v_rel ** 2)

    # values taken from Table 1 Potter & Collins (2013, MAPS).
    rho_i = 3314
    rho_t = 2629

    C_i = 5.43
    C_t = 3.816

    S_i = 1.34
    S_t = 1.28

    P_maxE = calc_Pmax(vimp_E, S_i, S_t, C_i, C_t, rho_i, rho_t) / 1e3
    P_maxL = calc_Pmax(vimp_L, S_i, S_t, C_i, C_t, rho_i, rho_t) / 1e3

    f_PE = 1 - (np.cos(0.5 * np.pi * 106 / (P_maxE * np.sin(45 * np.pi / 180)))) ** 1.3
    f_PL = 1 - (np.cos(0.5 * np.pi * 106 / (P_maxL * np.sin(45 * np.pi / 180)))) ** 1.3

    f_PE[np.where(P_maxE <= 106 * 2 / np.sqrt(2))] = 1
    f_PL[np.where(P_maxL <= 106 * 2 / np.sqrt(2))] = 1

    _ = plt.figure(figsize=(fig_width, fig_height))

    plt.plot(v_rel, f_PE, label='Earth', c=cm.bamako(0.3))
    plt.fill_between(v_rel, f_PE, alpha=0.25, color=cm.bamako(0.3))

    plt.plot(v_rel, f_PL, label='Moon', c=cm.bamako(0.6))
    plt.fill_between(v_rel, f_PL, f_PE, alpha=0.25, color=cm.bamako(0.6))

    # plt.text(2., 0.20, "unmelted")
    # plt.text(2., 0.85, "unmelted")

    plt.xlabel(r'$v_{\rm rel}$ [${\rm km}\,{\rm s}^{-1}$]', fontsize=13)
    plt.ylabel(r'Fraction of unmelted planetesimal', fontsize=13)

    plt.minorticks_on()

    plt.xlim(0, 40)
    plt.ylim(0, 1)

    plt.legend(fontsize=13, frameon=False)

    # plt.savefig('./figures/impactor_melting_vrel.pdf', format='pdf', bbox_inches='tight')

    plt.show()


def dunite_quartzite_comparison():
   
    v_rel = np.linspace(0, 50, 1000)

    vimp_E = np.sqrt(11.19 ** 2 + v_rel ** 2)
    vimp_L = np.sqrt(2.38 ** 2 + 3 + 1.022 ** 2 + v_rel ** 2)

    # values taken from Table 1 Potter & Collins (2013, MAPS).
    rho_i = 3314
    rho_t = 2629

    C_i = 5.43
    C_t = 3.816

    S_i = 1.34
    S_t = 1.28

    P_maxE = calc_Pmax(vimp_E, S_i, S_t, C_i, C_t, rho_i, rho_t) / 1e3

    f_PE_dun = 1 - (np.cos(0.5 * np.pi * 106 / (P_maxE * np.sin(45 * np.pi / 180)))) ** 1.3
    f_PE_dun[np.where(P_maxE <= 106 * 2 / np.sqrt(2))] = 1

    rho_i = 2649
    rho_t = 2629

    C_i = 1.143
    C_t = 3.816

    S_i = 1.43
    S_t = 1.28

    P_maxE = calc_Pmax(vimp_E, S_i, S_t, C_i, C_t, rho_i, rho_t) / 1e3

    f_PE_qua = 1 - (np.cos(0.5 * np.pi * 59 / (P_maxE * np.sin(45 * np.pi / 180)))) ** 1.3
    f_PE_qua[np.where(P_maxE <= 59 * 2 / np.sqrt(2))] = 1

    _ = plt.figure(figsize=(fig_width, fig_height))

    plt.plot(v_rel, f_PE_dun, label='Dunite', c=cm.bamako(0.6))
    plt.fill_between(v_rel, f_PE_qua, f_PE_dun, alpha=0.25, color=cm.bamako(0.6))


    plt.plot(v_rel, f_PE_qua, label='Quartzite', c=cm.bamako(0.3))
    plt.fill_between(v_rel, f_PE_qua, alpha=0.25, color=cm.bamako(0.3))

    plt.xlabel(r'$v_{\rm rel}$ [${\rm km}\,{\rm s}^{-1}$]', fontsize=13)
    plt.ylabel(r'Fraction of unmelted planetesimal', fontsize=13)

    plt.minorticks_on()

    plt.xlim(0, 40)
    plt.ylim(0, 1)

    plt.legend(fontsize=13, frameon=False)

    # plt.savefig('./figures/dunite_quartzite_comparison.pdf', format='pdf', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":

    plot_vimp_Pmax()
    main()
    dunite_quartzite_comparison()

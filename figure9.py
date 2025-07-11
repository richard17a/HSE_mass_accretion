import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import solve_ivp
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


def differential_equations(t: float, y: list, sigma_imp: float, rho_imp: float, eta: float, rho_atm0: float):
  """
  Defining the differential equations to be solved when calculating the impactors trajectory
  """

  ### ------ Defining constants ------
  C_d = 0.7 # the impactor's drag coefficient
  C_h = 0.02 # impactor's heat transfer coefficient (fraction of energy that heats impactor vs atmosphere)
  C_l = 0.001 # the impactor's lift coefficient
  M_E = 5.97e24
  R_E = 6371e3
  G = 6.67e-11
  sigma = 5.6704e-8
  T = 25000 # assumed temperatuer of shocked gas at leading edge of the impactor
  H = 8e3 # Earth's atmospheric scale height
  ### ----------------------- ###

  V, M, theta, Z, R, Rdot = y

  rho_a = rho_atm0 * np.exp(- Z / H)
  A = np.pi * R ** 2
  g = G * M_E / (R_E + Z) ** 2

  # equation 1 chyba et al. 1993 (I actually think there is a typo in their manuscript, and should read + g sin(theta))
  dVdt = - C_d * rho_a * A * V**2 / M + g * np.sin(theta)

  # equation 4 chyba et al. 1993
  dthetadt = (M * g * np.cos(theta) - 0.5 * C_l * rho_a * A * V**2) /\
              (M * V) - V * np.cos(theta) / (R_E + Z)

  # evaluate the evolution of impactor's altitude (vertical component of velocity vector)
  dZdt = - V * np.sin(theta)

  # equation 3 chyba et al. 1993
  dMdt = - np.minimum(sigma * T**4, 0.5 * C_h * rho_a * V**3) * A / eta

  if 0.25 * C_d * rho_a * V**2 > sigma_imp:
      R_dot = Rdot
      R_ddot = C_d * rho_a * V**2 / (2 * rho_imp * R)
  else:
      R_dot = 0.0
      R_ddot = 0.0

  return [dVdt, dMdt, dthetadt, dZdt, R_dot, R_ddot]


def event_Z_crossing(t: float, y: list):
  """
  Event triggered when altitude crosses 0 (i.e. the impactor hits the ground)
  """

  return y[3]


def event_mass_zero(t: float, y: list):
  """
  Event triggered when all mass has been ablated
  """

  return y[1]


def event_pancake(t: float, y: list, R0: float):
  """
  Event triggered when impactor's radius (size of pancake) exceeds 6 * initial radius
  At this point the model seems to break down - in reality individual fragments would
  form and develop their own bow shocks - so I just halt the integration at this point
  and classify it as an airburst event.
  """

  return 6 * R0 - y[4]


def event_dVdt_zero(t: float, y: list, rho_atm0: float):
  """
  Event triggered when object reaches terminal velocity
  """

  C_d = 0.7 # the impactor's drag coefficient
  H = 8e3 # Earth's atmospheric scale height
  G = 6.67e-11
  M_E = 5.97e24
  R_E = 6371e3

  V, M, theta, Z, R, _ = y
  rho_a = rho_atm0 * np.exp(- Z / H) # assuming an isothermal atmosphere
  A = np.pi * R ** 2 # surface area of leading edge of impactor
  g = G * M_E / (R_E + Z) ** 2

  term_vel = np.sqrt((M * g * np.sin(theta)) / (C_d * rho_a * A))

  out_num = 1

  if (np.isclose(term_vel, V, rtol=1e-01)) & (V < 1e3):
      """
      included the requirement for v<1km/s, as terminal velocity is initially
      very low due to low density in upper atmos
      """
      out_num = 0

  return out_num


def run_intergration(V0: float, M0: float, theta0: float, Z0: float, R0: float, Rdot0: float, sigma_imp: float, rho_imp: float, eta: float, rho_atm0=1.225):
  """
  This function will run the integration and calculate the atmospheric trjectory of the impactor

  Args:
      V0: Initial velocity of impactor
      M0: Initial mass of impactor
      theta0: Initial angle of impactor's trajectory
      Z0: Initial altitude of impactor
      R0: Intiial radius of impactor
      Rdot0: Initial rate of deformation (this is always zero...!)
      sigma_imp: Tensile strength of impactor
      rho_imp: Bulk density of impactor
      eta: impactor's heat of ablation
      rho_atmo0: Density of atmosphere at altitude=0
  """

  t_span = (0, 500)


  def event_pancake_with_R0(t: float, y: list):
      """
      Event triggered when size of pancake exceeds 6 * initial radius
      """

      return event_pancake(t, y, R0)


  def event_dVdt_zero_rhoatm0(t: float, y: list):
      """
      Event triggered when object reaches terminal velocity
      """

      return event_dVdt_zero(t, y, rho_atm0)

  event_Z_crossing.terminal = True
  event_Z_crossing.direction = -1

  event_mass_zero.terminal = True
  event_mass_zero.direction = -1

  event_dVdt_zero_rhoatm0.terminal = True
  event_dVdt_zero_rhoatm0.direction = 0

  event_pancake_with_R0.terminal = True
  event_pancake_with_R0.direction = -1

  # these events terminate the integration (i.e. when the impactor's mass = 0, or the altitude = 0 etc.)
  events = [event_Z_crossing, event_mass_zero, event_dVdt_zero_rhoatm0, event_pancake_with_R0]

  sol_iso = solve_ivp(
      fun=lambda t, y: differential_equations(t, y, sigma_imp, rho_imp, eta, rho_atm0),
      t_span=t_span,
      y0=[V0, M0, theta0, Z0, R0, Rdot0],
      method='RK45',
      dense_output=True,
      events=events,
      max_step=1e-2
  )

  # extract the solutions for the impactor's trajectory - i.e. all times until the integration is terminated
  t = sol_iso.t

  vel = sol_iso.sol(t)[0][:len(t)]
  mass = sol_iso.sol(t)[1][:len(t)]
  theta = sol_iso.sol(t)[2][:len(t)]
  altitude = sol_iso.sol(t)[3][:len(t)]
  radius = sol_iso.sol(t)[4][:len(t)]

  C_d = 0.7 # drag coefficient
  C_h = 0.02 # heat transfer coefficient
  M_E = 5.97e24
  R_E = 6371e3
  G = 6.67e-11
  sigma = 5.6704e-8
  T = 25000 # assumed temperature of shocked gas at leading edge of impactor
  H = 8e3 # Earth's atmospheric scale height

  g = G * M_E / (R_E + altitude) ** 2
  rho_a = rho_atm0 * np.exp(- altitude / H)
  A = np.pi * radius ** 2

  # equation 1 chyba et al. 1993 (I actually think there is a typo in their manuscript, and should read + g sin(theta))
  dVdt = - C_d * rho_a * A * vel**2 / mass + g * np.sin(theta)

  # equation 3 chyba et al. 1993
  dMdt = - np.minimum(sigma * T**4, 0.5 * C_h * rho_a * vel**3) * A / eta

  # evaluate the altitude evolution of the impactor's trajectory
  dZdt = - vel * np.sin(theta)

  ram_pressure = 0.7 * rho_atm0 * np.exp(-altitude / H) * vel ** 2 / 2

  # calculate the rate of change of the impactor's kinetic energy (wrt time, and altitude respectively)
  Ekindot = mass * vel * dVdt + 0.5 * vel**2 * dMdt
  dEkindh = Ekindot / dZdt

  return t, vel, mass, theta, altitude, radius, ram_pressure, dEkindh


def calc_Dmin(V0, rho_atm):

  rho_m = 3.5e3 # impactor's density
  Rdot0 = 0 # setting the initial rate of deformation of impactor to be 0 (required to solve 2nd order diff. eq. for radius expansion)
  theta0 = 45. * np.pi / 180. # initial angle of trajectory (most probable angle)
  Z0 = 100e3 # initial altitude of impctor

  R0 = 1
  M0 = rho_m * (4 * np.pi / 3) * (R0 ** 3) # initial mass of impactor

  melt_flag = False
  while melt_flag == False:

    _, velocity, _, _, altitude, _, _, _ = run_intergration(1.0 * V0, M0, theta0, Z0, R0, Rdot0, 1e6, rho_m, 8e6, rho_atm)

    if altitude[-1] > 10:
      R0 = 1.05 * R0
      M0 = rho_m * (4 * np.pi / 3) * (R0 ** 3)
    elif velocity[-1] < 15e3:
      R0 = 1.05 * R0
      M0 = rho_m * (4 * np.pi / 3) * (R0 ** 3)
    else:
      melt_flag = True

  return R0


def main():
   
    rho_m = 3.5e3 # impactor's density
    rho_atm = 1.225 # density of atmosphere at the surface (altitude = 0)
    V0 = 20e3 # initial velocity of impactor
    Rdot0 = 0 # this is just setting the initial rate of deformation of impactor to be 0 (required to solve 2nd order diff. eq. for radius expansion)
    theta0 = 45. * np.pi / 180. # initial angle of trajectory (most probable angle)
    Z0 = 100e3 # initial altitude of impactor  


    rho_atms = np.logspace(-2, 2, 8)

    D_mins  = []
    for rho_atm in rho_atms:

        D_mins = np.append(D_mins, calc_Dmin(20e3, rho_atm))

    _ = plt.figure(figsize=(fig_width, fig_height))

    plt.plot(rho_atms, D_mins, marker='.', c=cm.bamako(0.50), label=r'$v_{\rm entry} = 20\,{\rm km}\,{\rm s}^{-1}$')

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'$\rho_{\rm atm, \, surf.}\,[{\rm kg}\,{\rm m}^{-3}]$', fontsize=13)
    plt.ylabel('Critical impactor size for melting [m]', fontsize=13)

    plt.text(8, 350, r"$v_{\rm entry}=20\,{\rm km}\,{\rm s}^{-1}$", rotation=26, color=cm.bamako(0.25))
    plt.text(11, 175, r"$v_{\rm melt}=15\,{\rm km}\,{\rm s}^{-1}$", rotation=26, color=cm.bamako(0.25))

    plt.minorticks_on()

    plt.savefig('./figures/figure9.pdf', format='pdf', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
   
   main()

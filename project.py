#!/usr/local/bin/python

# Global constant, used to initialize the collider (GeV)
INITIAL_COLLIDER_ENERGY = 3

# Import 3rd party modules
import pylab as plb
from numpy import *
from math import ceil

# Import our modules
from lib.constants import *
from lib.functions import *
from lib.collider import Collider
  
#  Create a new collider and set its energy
collider = Collider(INITIAL_COLLIDER_ENERGY)

## BEGIN DIFFERENTIAL CROSS SECTIONS ##

"""
Example usage:
  gamma_gamma(0.1)
"""

def gamma_gamma(cos_theta):
  """Gamma-gamma differential cross section"""
  
  # This factor is used a couple of times, so capture it so as to
  # a) tidy things up and b) prevent wasted cycles
  epsilon_factor = 4*collider.epsilon_2
  
  alpha = (g_e**4 * sqrt(1 - epsilon_factor)) / (64*pi*pi*collider.s)
  
  # I'm using cos^2(theta) + sin^2(theta) = 1 here,
  # but this could be unstable if one term is small.
  # May only be unstable for recurrance relations?
  cos_theta_2 = cos_theta * cos_theta
  sin_theta_2 = 1 - cos_theta_2
  
  return alpha*(1 + cos_theta_2 + (epsilon_factor*sin_theta_2))
  
def z_z(cos_theta):
  """Z-Z differential cross section"""

  epsilon_factor = 4*collider.epsilon_2

  # Shorter name for neatness
  c = collider

  # Variables
  cos_theta_2 = cos_theta * cos_theta
  sin_theta_2 = 1 - cos_theta_2
  
  epsilon_factor = 1 - (4*c.epsilon_2)
  zeta_factor = 1 - (c.zeta*c.zeta)
  
  alpha   = (g_z**4 * sqrt(epsilon_factor)) / (1024*pi*pi*c.s)
  alpha  /= (zeta_factor*zeta_factor) + (c.zeta*c.zeta*gamma_z0*gamma_z0 / c.s)
  beta    = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_v)*(C_u_v)
  gamma   = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_a)*(C_u_a) * epsilon_factor
  delta   = 8 * C_e_v * C_e_a * C_u_v * C_u_a * sqrt(epsilon_factor)
  
  return alpha*(beta*(1 + cos_theta_2 + (4*c.epsilon_2*sin_theta_2)) + gamma*(1 + cos_theta_2) + delta*cos_theta)
  
def gamma_z(cos_theta):
  """gamma-Z differential cross section"""
  
  c = collider
  
  # Variables
  cos_theta_2 = cos_theta * cos_theta
  sin_theta_2 = 1 - cos_theta_2
  
  epsilon_factor = sqrt(1 - c.epsilon_2)
  zeta_factor = 1 - (c.zeta*c.zeta)
  
  alpha   = (2 * g_e*g_e * g_z*g_z * epsilon_factor * zeta_factor) / (256*pi*pi*c.s)
  alpha  /= (zeta_factor*zeta_factor) + (c.zeta*c.zeta * gamma_z0*gamma_z0  / c.s)
  beta    = C_e_v * C_u_v
  delta   = 2 * C_e_a * C_u_a * epsilon_factor
  
  return alpha * (beta*(1 + cos_theta_2 + (4*c.epsilon_2*sin_theta_2)) + (delta*cos_theta))
  
def combined_diff_cross_section(cos_theta):
  """Sum of gamma-gamma, z-z and gamma-z differential cross sections."""
  return gamma_gamma(cos_theta) + z_z(cos_theta) + gamma_z(cos_theta)
  
## END DIFFERENTIAL CROSS SECTIONS ##


## BEGIN DIFFERENTIAL PLOTS ##

"""
Example usage:
  plot_gamma_gamma(-1, 1)
"""

def plot_gamma_gamma(a, b, step_size = 0.1):
  """Plot gamma-gamma differential cross section between a and b"""
  
  cos_theta = arange(a, b+step_size, step_size)
  z = [gamma_gamma(x) for x in cos_theta]
  plb.plot(cos_theta, z)
  plb.show()


def plot_z_z(a, b, step_size = 0.1):
  """Plot z-z differential cross section between a and b"""

  cos_theta = arange(a, b+step_size, step_size)
  z = [z_z(x) for x in cos_theta]
  plb.plot(cos_theta, z)
  plb.show()

def plot_gamma_z(a, b, step_size = 0.1):
  """Plot gamma-z differential cross section between a and b"""

  cos_theta = arange(a, b+step_size, step_size)
  z = [gamma_z(x) for x in cos_theta]
  print [cos_theta[0], z[0]], [cos_theta[-1], z[-1]]
  plb.plot(cos_theta, z)
  plb.show()
 
def plot_combined(a, b, step_size = 0.1):
  """Plot the sum of all the differential cross sections between a and b"""
  
  cos_theta = arange(a, b, step_size)
  z = [combined_diff_cross_section(x) for x in cos_theta]
  plb.plot(cos_theta, z)
  plb.show()

## END DIFFERENTIAL PLOTS ##


## BEGIN CONSISTENCY CHECKS ##

"""
Example usage:
  check_gamma_gamma_consistency()
"""

def check_gamma_gamma_consistency():
  """Checks to see if the numerical gamma-gamma cross section complies with theory"""

  # Sample plot from cos(theta) = -1..1
  # plot_gamma_gamma(-1, 1)

  # Expected and calculated values
  factor = 4*collider.epsilon_2
  alpha = (g_e**4 / (64*pi*pi*collider.s)) * sqrt(1 - factor) 
  print "\nExpected value:    {0:.6}".format((4.0/3.0)*alpha*(2 + factor))
  print "Trapezium value:   {0:.6}".format(trapezium(gamma_gamma, -1, 1, 10000))
  print "Monte carlo value: {0:.6}".format(montecarlo2(gamma_gamma, -1, 1, 1000))

  # N = 100:
  # Expected value:  4.52108e-06
  # Trapezium value: 4.52131e-06
  # N = 1000:
  # Expected value:  4.52108e-06
  # Trapezium value: 4.52108e-06
  
def check_z_z_consistency():
  """Checks to see if the numerical z-z cross section complies with theory"""

  # Sample plot from cos(theta) = -1..1
  # plot_z_z(-1, 1)

  # Expected and calculated values
  # See notes for greek letter translations
  c = collider
  
  epsilon_factor = 1 - (4*c.epsilon_2)
  zeta_factor = 1 - (c.zeta*c.zeta)
  
  alpha   = (g_z**4 * sqrt(epsilon_factor)) / (1024*pi*pi*c.s)
  alpha  /= (zeta_factor*zeta_factor) + (c.zeta*c.zeta*gamma_z0*gamma_z0 / c.s)
  beta    = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_v)*(C_u_v)
  gamma   = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_a)*(C_u_a) * epsilon_factor
  # No need for delta as it isn't a term in sigma

  # Haven't integrated over phi, so no 2*pi prefactor
  print "\nExpected value:    {0:.6}".format((8.0/3.0)*alpha*(gamma + ((1 + 2*c.epsilon_2)*beta)))
  print "Trapezium value:   {0:.6}".format(trapezium(z_z, -1, 1, 1000))
  print "Monte carlo value: {0:.6}".format(montecarlo2(z_z, -1, 1, 10000))

  # N = 100:
  # Expected value:    6.56796e-13
  # Trapezium value:   6.56828e-13
  # N = 1000:
  # Expected value:    6.56796e-13
  # Trapezium value:   6.56796e-13  

def check_gamma_z_consistency():
  """Checks to see if the numerical gamma-z cross section complies with theory"""
  
  # Sample plot from cos(theta) = -1..1
  # plot_gamma_z(-1, 1)

  # Expected and calculated values
  # See notes for greek letter translations
  c = collider
  epsilon_factor = sqrt(1 - c.epsilon_2)
  zeta_factor = 1 - (c.zeta*c.zeta)
  
  alpha   = (2 * g_e*g_e * g_z*g_z * epsilon_factor * zeta_factor) / (256*pi*pi*c.s)
  alpha  /= (zeta_factor*zeta_factor) + (c.zeta*c.zeta * gamma_z0*gamma_z0  / c.s)
  beta    = C_e_v * C_u_v
  delta   = 2 * C_e_a * C_u_a * epsilon_factor

  print "\nExpected value:    {0:.6}".format((8.0/3.0)*alpha*beta*(1 + 2*c.epsilon_2))
  print "Trapezium value:   {0:.6}".format(trapezium(gamma_z, -1, 1, 1000))
  # mc only starts to get reasonable results at 10e-6 points
  print "Monte carlo value: {0:.6}".format(montecarlo2(gamma_z, -1, 1, 1000))

## END CONSISTENCY CHECKS ##


## BEGIN THEORY CROSS SECTIONS ##

"""
Example usage:
  theory_root_s, theory_cross_section  = gamma_gamma_theory_cross_section(3, 30)
  plb.plot(theory_root_s, theory_cross_section)
  plb.show()
"""

def gamma_gamma_theory_cross_section(a, b, step_size = 0.1):
  """Evaluate between collider energy a to b in steps of step_size"""
  
  # A range of collider energies from 3GeV to 30GeV
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  # Array to store the integrated differential cross section, sigma
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    
    factor = 4*collider.epsilon_2
    alpha = (g_e**4 / (64*pi*pi*collider.s)) * sqrt(1 - factor) 
    
    # The theoretical cross section from integrating the (simple) diff. cross section
    cross_section[count] = (4.0/3.0)*alpha*(2 + factor)
    count += 1

  return root_s_arr, cross_section

def z_z_theory_cross_section(a, b, step_size = 0.1):
  """Evaluate between collider energy a to b in steps of step_size"""
  # This one's is good to go from 3 to 200GeV. Yeeeah.
  
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    c = collider

    epsilon_factor = 1 - (4*c.epsilon_2)
    zeta_factor = 1 - (c.zeta*c.zeta)

    alpha   = (g_z**4 * sqrt(epsilon_factor)) / (1024*pi*pi*c.s)
    alpha  /= (zeta_factor*zeta_factor) + (c.zeta*c.zeta*gamma_z0*gamma_z0 / c.s)
    beta    = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_v)*(C_u_v)
    gamma   = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_a)*(C_u_a) * epsilon_factor

    cross_section[count] = (8.0/3.0)*alpha*(gamma + ((1 + 2*c.epsilon_2)*beta))
    count += 1
    
  return root_s_arr, cross_section
  
def gamma_z_theory_cross_section(a, b, step_size = 0.1):
  """Evaluate between collider energy a to b in steps of step_size"""
  
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    c = collider
    
    epsilon_factor = sqrt(1 - c.epsilon_2)
    zeta_factor = 1 - (c.zeta*c.zeta)

    alpha   = (2 * g_e*g_e * g_z*g_z * epsilon_factor * zeta_factor) / (256*pi*pi*c.s)
    alpha  /= (zeta_factor*zeta_factor) + (c.zeta*c.zeta * gamma_z0*gamma_z0  / c.s)
    beta    = C_e_v * C_u_v
    delta   = 2 * C_e_a * C_u_a * epsilon_factor
    
    cross_section[count] = (8.0/3.0) * alpha * beta * (1 + 2*c.epsilon_2)
    count += 1
    
  return root_s_arr, cross_section
  
## END THEORY CROSS SECTIONS ##

## BEGIN NUMERICAL CROSS SECTIONS ##

"""
Example usage:
  mc_root_s, mc_cross_section  = montecarlo_cross_section(gamma_gamma, 3, 30)
  plb.plot(mc_root_s, mc_cross_section)
  plb.show()
"""

def trapezium_cross_section(f, a, b, step_size = 0.1, strips = 1000):
  """Compete the cross section using the trapezium rule in a range of collider energies."""
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    # Use the trapezium rule to calculate the numerical cross section
    cross_section[count] = trapezium(f, -1, 1, strips)
    count += 1

  return root_s_arr, cross_section


def montecarlo_cross_section(f, a, b, step_size = 0.1, iterations = 1000):
  """Compete the cross section using monte carlo integration in a range of collider energies."""
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    # Use the monte carlo method to calculate the numerical cross section
    cross_section[count] = montecarlo2(f, -1, 1, iterations)
    count += 1
    
  return root_s_arr, cross_section

## END NUMERICAL CROSS SECTIONS ##


## BEGIN CROSS SECTION COMPARISONS ##

"""
Example usage:
  compare_gamma_gamma(3, 30)
"""

def compare_gamma_gamma(a, b, step_size = 0.1):
  """Compare the theoretical with the numerical cross sections."""
  theory_root_s, theory_sigma = gamma_gamma_theory_cross_section(a, b, step_size)
  trap_root_s, trap_sigma = trapezium_cross_section(gamma_gamma, a, b, step_size)
  mc_root_s, mc_sigma = montecarlo_cross_section(gamma_gamma, a, b, step_size)

  plb.plot(theory_root_s, theory_sigma)
  plb.plot(trap_root_s, trap_sigma)
  plb.plot(mc_root_s, mc_sigma)
  # Pretty good fit tbh
  plb.show()
  
def compare_z_z(a, b, step_size = 0.1):
  theory_root_s, theory_sigma = z_z_theory_cross_section(a, b, step_size)
  trap_root_s, trap_sigma = trapezium_cross_section(z_z, a, b, step_size)
  mc_root_s, mc_sigma = montecarlo_cross_section(z_z, a, b, step_size)

  plb.plot(theory_root_s, theory_sigma)
  plb.plot(trap_root_s, trap_sigma)
  plb.plot(mc_root_s, mc_sigma)
  plb.show()
  
def compare_gamma_z(a, b, step_size = 0.1):
  theory_root_s, theory_sigma = gamma_z_theory_cross_section(a, b, step_size)
  trap_root_s, trap_sigma = trapezium_cross_section(gamma_z, a, b, step_size)
  mc_root_s, mc_sigma = montecarlo_cross_section(gamma_z, a, b, step_size, 1000)

  plb.plot(theory_root_s, theory_sigma)
  plb.plot(trap_root_s, trap_sigma)
  plb.plot(mc_root_s, mc_sigma)
  plb.show()
  
def compare_gamma_gamma_to_z_z(a, b, step_size = 0.1):
  """Compare the theoretical gamma-gamma and z-z cross sections."""
  # g_g_root_s, g_g_cross_section  = gamma_gamma_theory_cross_section(a, b, step_size)
  # z_z_root_s, z_z_cross_section  = z_z_theory_cross_section(a, b, step_size)
  g_g_root_s, g_g_cross_section  = montecarlo_cross_section(gamma_gamma, a, b, step_size)
  z_z_root_s, z_z_cross_section  = montecarlo_cross_section(z_z, a, b, step_size)

  combined_theory_cross_section = []
  for i in range(len(z_z_cross_section)):
    combined_theory_cross_section.append(g_g_cross_section[i] + z_z_cross_section[i])

  # Either root_s array will work here, same scale
  plb.plot(g_g_root_s, combined_theory_cross_section)
  # 3-300GeV, jesus christ that is beautiful
  plb.show()

def compare_all(a, b, step_size = 0.1):
  """Compare all the theoretical (g-g, z-z, g-z) crosss sections."""
  g_g_root_s, g_g_cross_section  = gamma_gamma_theory_cross_section(a, b, step_size)
  z_z_root_s, z_z_cross_section  = z_z_theory_cross_section(a, b, step_size)
  g_z_root_s, g_z_cross_section  = gamma_z_theory_cross_section(a, b, step_size)
  
  combined_theory_cross_section_1 = []
  combined_theory_cross_section_2 = []
  for i in range(len(z_z_cross_section)):
    # Cross section with, then without, the interference term
    combined_theory_cross_section_1.append(g_g_cross_section[i] + z_z_cross_section[i])
    combined_theory_cross_section_2.append(g_g_cross_section[i] + z_z_cross_section[i] + g_z_cross_section[i])
    
  plb.plot(g_g_root_s, combined_theory_cross_section_1)
  plb.plot(g_g_root_s, combined_theory_cross_section_2)
  plb.show()
  
## END CROSS SECTION COMPARISONS ##


## BEGIN NUMERICAL CROSS SECTION PLOTS ##

"""Example usage:
  root_s_range, sigma = seperate_cross_section('trapezium', 3, 100)
  plb.plot(root_s_range, sigma)
  plb.show()
"""

method_map = {
  'trapezium':  trapezium_cross_section,
  'mc':         montecarlo_cross_section,
  'montecarlo': montecarlo_cross_section
}

def seperate_cross_section(method, a, b, step_size = 0.1, N = 1000):
  """
  Integrate each cross section seperately in a range of collider energies, then add the result.
  N is number of strips for trapezium, number of points for monte carlo.
  """
  f = method_map[method]
  g_g_root_s, g_g_sigma = f(gamma_gamma, a, b, step_size, N)
  z_z_root_s, z_z_sigma = f(z_z, a, b, step_size, N)
  g_z_root_s, g_z_sigma = f(gamma_z, a, b, step_size, N)
  
  # Could use the length of any list here, they all have the same dimension
  length = len(g_g_root_s)
  combined_cross_section = zeros(length)
  for i in range(length):
    combined_cross_section[i] = g_g_sigma[i] + z_z_sigma[i] + g_z_sigma[i]
  
  return g_g_root_s, combined_cross_section
  
def combined_cross_section(method, a, b, step_size = 0.1, N = 1000):
  """Integrate the combined cross section in a range of collider energies."""
  f = method_map[method]
  return f(combined_diff_cross_section, a, b, step_size, N)
  
## END NUMERICAL CROSS SECTION PLOTS ##
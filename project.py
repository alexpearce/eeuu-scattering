#!/usr/local/bin/python

# Global constant, used to initialize the collider
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

## BEGIN GAMMA-GAMMA ##

def gamma_gamma(cos_theta):
  """Gamma-gamma differential cross section"""
  
  # This factor is used a couple of times, so capture it so as to
  # a) tidy things up and b) prevent wasted cycles
  epsilon_factor = 4*collider.epsilon_2
  
  coefficient = (g_e**4 / ((8*pi)*(8*pi)*collider.s)) * sqrt(1 - epsilon_factor)
  
  # I'm using cos^2(theta) + sin^2(theta) = 1 here,
  # but this could be unstable if one term is small.
  # May only be unstable for recurrance relations?
  cos_theta_2 = cos_theta * cos_theta
  sin_theta_2 = 1 - cos_theta_2
  
  variable = 1 + cos_theta_2 + (epsilon_factor*sin_theta_2)
  
  return coefficient*variable
  
def plot_gamma_gamma(a, b, step_size = 0.1):
  """Plot gamma-gamma differential cross section between a and b"""
  
  cos_theta = arange(a, b+step_size, step_size)
  z = [gamma_gamma(x) for x in cos_theta]
  plb.plot(cos_theta, z)
  plb.show()
  
def check_gamma_gamma_consistency():
  """Checks to see if the numerical gamma-gamma cross section complies with theory"""

  # Sample plot from cos(theta) = -1..1
  # plot_gamma_gamma(-1, 1)

  # Expected and calculated values
  factor = 4*collider.epsilon_2
  coefficient = (g_e**4 / ((8*pi)*(8*pi)*collider.s)) * sqrt(1 - factor) 
  print "Expected value:    {0:.6}".format(2*coefficient*((1 + factor) +(1.0/3.0)*(1 - factor)))
  print "Trapezium value:   {0:.6}".format(trapezium(gamma_gamma, -1, 1, 1000))
  print "Monte carlo value: {0:.6}".format(montecarlo(gamma_gamma, 10000))

  # N = 100:
  # Expected value:  4.52108e-06
  # Trapezium value: 4.52131e-06
  # N = 1000:
  # Expected value:  4.52108e-06
  # Trapezium value: 4.52108e-06

# check_gamma_gamma_consistency()
  
## END GAMMA-GAMMA ##

## BEGIN Z-Z ##
  
def z_z(cos_theta):
  """Z-Z differential cross section"""

  epsilon_factor = 4*collider.epsilon_2
  
  # Shorter name for neatness
  c = collider
  
  # Variables
  cos_theta_2 = cos_theta * cos_theta
  sin_theta_2 = 1 - cos_theta_2
  
  # We'll do this line by line, similarly to the project description
  # Note that it is *not* the case that all lines are simply multiplied together
  line_1_numerator = (g_z**4 * sqrt(1 - epsilon_factor)) / (1024*pi*pi*collider.s)
  line_1_denominator = (1 - 2*c.lamba*c.lamba + c.lamba**4) + ((c.lamba*gamma_z)*(c.lamba*gamma_z) / c.s)
  line_1 = line_1_numerator / line_1_denominator
  
  line_2 = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_v)*(C_u_v) * (1 + cos_theta_2 + (epsilon_factor*sin_theta_2))
  line_3 = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_a)*(C_u_a) * (1 + cos_theta_2) * (1 - epsilon_factor)
  line_4 = 8 * C_e_v * C_e_a * C_u_v * C_u_a * sqrt(1 - epsilon_factor) * cos_theta
  
  return line_1 * (line_2 + line_3 + line_4)
  
def plot_z_z(a, b, step_size = 0.1):
  """Plot z-z differential cross section between a and b"""

  cos_theta = arange(a, b+step_size, step_size)
  z = [z_z(x) for x in cos_theta]
  plb.plot(cos_theta, z)
  plb.show()

def check_z_z_consistency():
  """Checks to see if the numerical z-z cross section complies with theory"""

  # Sample plot from cos(theta) = -1..1
  # plot_z_z(-1, 1)

  # Expected and calculated values
  # See written notes for greek letter translations
  # TODO: latex these
  c = collider
  epsilon_plus = 1 + (4*c.epsilon*c.epsilon)
  epsilon_minus = 1 - (4*c.epsilon*c.epsilon)
  alpha   = (g_z**4 * sqrt(epsilon_minus)) / (1024*pi*pi*c.s)
  alpha  /= (1 - (c.lamba*c.lamba))**2 + (c.lamba*c.lamba*gamma_z*gamma_z / c.s)
  beta    = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_v)*(C_u_v)
  delta   = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_a)*(C_u_a)
  
  # Haven't integrated over phi, no no 2*pi prefactor
  prefactor   = 2*alpha
  beta_term   = (epsilon_plus + (1.0/3.0)*epsilon_minus) * beta
  delta_term  = (4.0/3.0) * epsilon_minus * delta
  
  print "Expected value:    {0:.6}".format(prefactor * (beta_term + delta_term))
  print "Trapezium value:   {0:.6}".format(trapezium(z_z, -1, 1, 1000))
  # TODO: fix monte carlo integration (it's hardcoded on gamma-gamma).
  # print "Monte carlo value: {0:.6}".format(montecarlo(z_z, 1000))

  # N = 100:
  # Expected value:    6.56796e-13
  # Trapezium value:   6.56828e-13
  # N = 1000:
  # Expected value:    6.56796e-13
  # Trapezium value:   6.56796e-13

check_z_z_consistency()

## END Z-Z ##

## BEGIN GAMMA-Z ##
  
def gamma_z(cos_theta):
  """gamma-Z differential cross section"""
  pass

## END GAMMA-Z ##

def theory_cross_section(a, b, step_size = 0.1):
  """Evaluate between collider energy a to b in steps of step_size"""
  
  # A range of collider energies from 3GeV to 30GeV
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  # Array to store the integrated differential cross section, sigma
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    
    factor = 4*collider.epsilon_2
    coefficient = (g_e**4 / ((8*pi)*(8*pi)*collider.s)) * sqrt(1 - factor)
    # The theoretical cross section from integrating the (simple) diff. cross section
    cross_section[count] = 2*coefficient*((1 + factor) +(1.0/3.0)*(1 - factor))
    count += 1

  return root_s_arr, cross_section

def trapezium_cross_section(a, b, step_size = 0.1):
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    # Use the trapezium rule to calculate the numerical cross section
    cross_section[count] = trapezium(gamma_gamma, -1, 1, 1000)
    count += 1

  return root_s_arr, cross_section
  
def montecarlo_cross_section(a, b,step_size = 0.1):
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    collider.set_energy_to(i)
    # Use the monte carlo method to calculate the numerical cross section
    cross_section[count] = montecarlo(gamma_gamma, 1000)
    count += 1
    
  return root_s_arr, cross_section
  
# mc_root_s, mc_sigma = montecarlo_cross_section(3, 30)
# plb.plot(mc_root_s, mc_sigma)
# plb.show()

def compare_gamma_gamma(a, b, step_size = 0.1):
  theory_root_s, theory_sigma = theory_cross_section(a, b, step_size)
  trap_root_s, trap_sigma = trapezium_cross_section(a, b, step_size)
  mc_root_s, mc_sigma = montecarlo_cross_section(a, b, step_size)

  plb.plot(theory_root_s, theory_sigma)
  plb.plot(trap_root_s, trap_sigma)
  plb.plot(mc_root_s, mc_sigma)
  # Pretty good fit tbh
  plb.show()
  
# compare_gamma_gamma(3, 30)
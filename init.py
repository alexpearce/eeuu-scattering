#!/usr/local/bin/python

INITIAL_COLLIDER_ENERGY = 3

import pylab as plb
from numpy import *
from math import ceil

from lib.constants import *
from lib.functions import *
# Import my collider..
from lib.collider import Collider
  
#  Create a new collider and set its energy
c = Collider(INITIAL_COLLIDER_ENERGY)

# Gamma-gamma differential cross section
def gamma_gamma(cos_theta):
  # Uses the *global* root_s value
  factor = 4*c.epsilon_2
  
  coefficient = (g_e**4 / ((8*pi)*(8*pi)*c.s)) * sqrt(1 - factor)
  
  # I'm using cos^2(theta) + sin^2(theta) = 1 here,
  # but this could be unstable if one term is small.
  # May only be unstable for recurrance relations?
  cos_theta_2 = cos_theta * cos_theta
  sin_theta_2 = 1 - cos_theta_2
  
  variable = 1 + cos_theta_2 + (factor*sin_theta_2)
  
  return coefficient*variable
  
# Plot gamma-gamma differential cross section between a and b
def plot_gamma_gamma(a, b, step_size = 0.1):
  cos_theta = arange(a, b+step_size, step_size)
  z = [gamma_gamma(x) for x in cos_theta]
  plb.plot(cos_theta, z)
  plb.show()
  
# Z-Z differential cross section
def z_z(cos_theta):
  pass
  
# gamma-Z differential cross section
def gamma_z(cos_theta):
  pass

# Checks to see if the numerical gamma-gamma cross section complies with theory
def check_consistency():
  # Sample plot from cos(theta) = -1..1
  # plot_gamma_gamma(-1, 1)

  # Expected and calculated values
  factor = 4*c.epsilon_2
  coefficient = (g_e**4 / ((8*pi)*(8*pi)*c.s)) * sqrt(1 - factor) 
  print "Expected value:    {0:.6}".format(2*coefficient*((1 + factor) +(1.0/3.0)*(1 - factor)))
  print "Trapezium value:   {0:.6}".format(trapezium(gamma_gamma, -1, 1, 1000))
  print "Monte carlo value: {0:.6}".format(montecarlo(gamma_gamma, 10000))

  # N = 100:
  # Expected value:  4.52108e-06
  # Trapezium value: 4.52131e-06
  # N = 1000:
  # Expected value:  4.52108e-06
  # Trapezium value: 4.52108e-06
  
# check_consistency()

# Evaluate between collider energy a to b in steps of step_size
def theory_cross_section(a, b, step_size = 0.1):
  # A range of collider energies from 3GeV to 30GeV
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  # Array to store the integrated differential cross section, sigma
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    c.set_energy_to(i)
    
    factor = 4*c.epsilon_2
    coefficient = (g_e**4 / ((8*pi)*(8*pi)*c.s)) * sqrt(1 - factor)
    # The theoretical cross section from integrating the (simple) diff. cross section
    cross_section[count] = 2*coefficient*((1 + factor) +(1.0/3.0)*(1 - factor))
    count += 1

  return root_s_arr, cross_section

def trapezium_cross_section(a, b, step_size = 0.1):
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    c.set_energy_to(i)
    # Use the trapezium rule to calculate the numerical cross section
    cross_section[count] = trapezium(gamma_gamma, -1, 1, 1000)
    count += 1

  return root_s_arr, cross_section
  
def montecarlo_cross_section(a, b,step_size = 0.1):
  root_s_arr = [i for i in arange(float(a), float(b), step_size)]
  cross_section = zeros(len(root_s_arr))
  count = 0
  for i in root_s_arr:
    c.set_energy_to(i)
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
  
compare_gamma_gamma(3, 30)
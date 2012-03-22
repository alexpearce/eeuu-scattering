from math import *

# Particle data found at the PDG
# http://pdg.lbl.gov/2011/tables/contents_tables.html

# All units given in GeV, c = 1

### Constants ###

# Fine structure
alpha = 1.0 / 128.0

# Weak mixing angle / Weinberg angle
sin_2_theta_w = 0.23152
theta_w = asin(sqrt(sin_2_theta_w))

# QED e- coupling
g_e = sqrt(4*pi*alpha)

# QED Z coupling
g_z = g_e / (cos(theta_w) * sin(theta_w))

# Muon rest mass
m_u = 0.105658

# Z0 rest mass
m_z = 91.1876

# Z0 width
gamma_z0 = 2.5

# e- vectorial coupling
C_e_v = -0.5 + (2*sin_2_theta_w)
# e- axial coupling
C_e_a = -0.5

# mu vectorial coupling
C_u_v = C_e_v
# mu axial coupling
C_u_a = C_e_a

### End Constants ###
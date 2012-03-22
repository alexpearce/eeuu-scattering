from lib.constants import *

class Collider:
  """A few functions which calculate collider-energy-dependant quantities."""
  
  def __init__(self, energy):
    """`Collider(5)` creates a 5GeV collider"""
    self.set_energy_to(energy)
    
  def __repr__(self):
    """`print collider_instance` will show this string."""
    return "Current collider energy sqrt(s): {:.2} GeV".format(float(self.energy))
    
  def set_energy_to(self, energy):
    """Update all s-dependant quantities"""
    # Collider energy is sqrt(s) in GeV
    self.energy = self.root_s = float(energy)
    self.s      = self.root_s * self.root_s
    self.set_epsilon(self.energy)
    self.set_lamba(self.energy)
    
  def set_epsilon(self, energy):
    """Dimensionless variable. In practice we only use the square of epsilon."""
    self.epsilon = m_u / energy
    self._set_epsilon_2(self.epsilon)
    
  def _set_epsilon_2(self, epsilon):
    self.epsilon_2 = epsilon * epsilon
    
  def set_lamba(self, energy):
    """Dimensionless variable."""
    self.lamba = m_z / energy
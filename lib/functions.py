from numpy import random
from lib.constants import *
  
def trapezium(f, a, b, N):
  """
  Trapezium integration of a function f between a and b
  with N strips.
  """
  
  # Width of strips
  h = (b - a) / float(N)
  
  prefactor = h/2.0 * (f(a) + f(b))
  
  # Sum
  s = 0
  for i in range(1, N):
    s += f(a + i*h)
  
  return prefactor + h*s
  
def montecarlo(f, N):
  """
  Monte carlo integration of a function f between -1 and 1
  with N trials.
  """
  
  # 1. Generate a random (x, y) pair
  # 2. Compute f(x)
  # 3. If y < f(x), (x, y) is under the curve
  # 4. Compute the ratio hits / N
  # 5. Multiply the ratio by the box area
  hits = 0
  for i in range(N-1):
    # Generate x between -1 and 1
    x = 2*random.random() - 1
    # Generate y between 0 and 0.0000034
    y = 3.5e-6*random.random()
    if y < f(x): hits += 1
  
  hits_ratio = float(hits) / float(N)
  integration_area = 2.0 * 3.5e-6

  return integration_area * hits_ratio
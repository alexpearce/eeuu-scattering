from numpy import random, arange
import csv
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
  
def montecarlo(f, a, b, N):
  """
  Monte carlo integration of a function f between a and b (a < b)
  with N trials.
  No longer in use, use montecarlo2 instead.
  """
  
  # 1. Generate a random (x, y) pair
  # 2. Compute f(x)
  # 3. If y < f(x), (x, y) is under the curve
  # 4. Compute the ratio hits / N
  # 5. Multiply the ratio by the box area
  # Hits above and below the x-axis
  positive_hits = 0
  negative_hits = 0
  # The maximum and minimum value f takes between -1 and 1
  maximum, minimum = max_min_of(f, -1, 1)
  for i in range(N):
    # Generate x between -1 and 1
    x = (b - a)*random.random_sample() + a
    # Generate y between 0 and the function's maximum value
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.random.html#numpy.random.random
    y = (maximum - minimum)*random.random_sample() + minimum
    
    # Store the function value at x
    f_x = f(x)
    
    # y needs to be between the x-axis and the curve
    if f_x > 0:
      if 0 <= y <= f_x: positive_hits += 1
    elif f_x < 0:
      if f_x <= y <= 0: negative_hits += 1
    # else f(x) = 0, no area to integrate.
  
  hits_ratio = float(positive_hits - negative_hits) / float(N)
  # Integrating in a box of width 2 (between -1 and 1), height maximum - minimum
  integration_area = (b - a) * (maximum - minimum)

  return integration_area * hits_ratio
  
def montecarlo2(f, a, b, N):
  """A considerably better mc method. Same API as the first."""
  prefactor = float(b - a) / float(N)
  thesum = 0
  
  for i in range(N):
    x = (b - a)*random.random_sample() + a
    thesum += f(x)
  
  return prefactor * thesum
  
def max_min_of(f, a, b, step_size = 0.01):
  """Find the maximum value of a function f between a and b"""
  maximum = minimum = f(a)
  for x in arange(a, b+step_size, step_size):
    if (f(x) > maximum): maximum = f(x)
    if (f(x) < minimum): minimum = f(x)
  
  return maximum, (minimum if minimum < 0 else 0) 
  
def export_to_csv(x, y, filename = 'data'):
  """Exports two equal-dimension arrays to a CSV file"""
  if (len(x) != len(y)):
    print "x and y must be of the same dimension."
    exit()
  print "Writing 2D data to file..."
  writer = csv.writer(open('{}.csv'.format(filename), 'wb'), delimiter = ',')
  for k, v in enumerate(x):
    writer.writerow([v, y[k]])

  print "File writing complete."
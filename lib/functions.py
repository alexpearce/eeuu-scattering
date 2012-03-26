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
  # Hits above and below the x-axis
  positive_hits = 0
  negative_hits = 0
  # The maximum and minimum value f takes between -1 and 1
  maximum, minimum = max_min_of(f, -1, 1)
  x_arr = []
  y_arr = []
  print maximum, minimum
  for i in range(N-1):
    # Generate x between -1 and 1
    x = 2*random.random() - 1
    x_arr.append(x)
    # Generate y between 0 and the function's maximum value
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.random.html#numpy.random.random
    y = (maximum - minimum)*random.random() + minimum
    y_arr.append(y)
    if f(x) < 0:
      if f(x) < y < 0: negative_hits += 1
    elif f(x) > 0:
      if 0 < y < f(x): positive_hits += 1
    # else f(x) = 0, no area to integrate.
  
  print positive_hits, negative_hits
  positive_hits_ratio = float(positive_hits) / float(N)
  negative_hits_ratio = float(negative_hits) / float(N)
  # Integrating in a box of width 2 (between -1 and 1), height maximum - minimum
  integration_area = 2.0 * (maximum - minimum)

  return integration_area * (positive_hits_ratio - negative_hits_ratio)
  
def max_min_of(f, a, b, step_size = 0.001):
  """Find the maximum value of a function f between a and b"""
  maximum = minimum = f(a)
  for x in arange(a, b+step_size, step_size):
    if (f(x) > maximum): maximum = f(x)
    if (f(x) < minimum): minimum = f(x)
  
  return maximum, (minimum if minimum < 0 else 0) 

# TODO: export plots as PDF
def export_to_pdf():
  """docstring for export_as_data"""
  pass
  
def export_to_csv(x, y, filename = 'data'):
  """Exports two equal-dimension arrays to a CSV file"""
  if (len(x) != len(y)):
    print "x and y must be of the same dimension."
    exit()
  print "Writing data to file..."
  writer = csv.writer(open('{}.csv'.format(filename), 'wb'), delimiter = ' ')
  writer.writerow(['x', 'y'])
  for k, v in enumerate(x):
    writer.writerow([v, y[k]])

  print "File writing complete."
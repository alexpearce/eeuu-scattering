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
  hits = 0
  # The maximum value f takes between -1 and 1
  maximum = max_of(f, -1, 1)
  for i in range(N-1):
    # Generate x between -1 and 1
    x = 2*random.random() - 1
    # Generate y between 0 and the function's maximum value
    y = maximum*random.random()
    if y < f(x): hits += 1

  hits_ratio = float(hits) / float(N)
  # Integrating in a box of width 2 (between -1 and 1), height maximum
  integration_area = 2.0 * maximum

  return integration_area * hits_ratio
  
def max_of(f, a, b, step_size = 0.1):
  """Find the maximum value of a function f between a and b"""
  maximum = f(a)
  for x in arange(a+step_size, b+step_size, step_size):
    if (f(x) > maximum): maximum = f(x)
  
  return maximum
  
# TODO: export as PDF, export as CSV/some data format

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
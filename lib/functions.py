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
  s = sum([f(a + i*h) for i in range(1, N)])
  
  return prefactor + h*s
  
def montecarlo(f, a, b, N):
  """A considerably better mc method. Same API as the first."""
  prefactor = float(b - a) / float(N)
  thesum = 0
  
  difference = b - a
  for i in range(N):
    # Generate random number in the range [a,b), a > b.
    x = difference*random.random_sample() + a
    thesum += f(x)
  
  return prefactor * thesum
  
def export_to_csv(x, y, filename = 'data'):
  """Exports two equal-dimension arrays to a CSV file"""
  if (len(x) != len(y)):
    print "x and y must be of the same dimension."
    exit()
  print "Writing 2D data to {}.csv...".format(filename)
  writer = csv.writer(open('{}.csv'.format(filename), 'wb'), delimiter = ',')
  for k, v in enumerate(x):
    writer.writerow([v, y[k]])

  print "File writing complete."
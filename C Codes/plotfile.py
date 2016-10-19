"""Plotfile.py: small utility to support SESG6025.

Usage:

python plotfile.py FILENAME

expects data in the ascii data file FILENAME, organised in row and columns, and
attempts to plot column 2 against column 1 if there are two columns, or
attempts to plot columns 2 to n against column 1 if there are n columns (n>2)

Hans Fangohr, Oct 2010, revised Oct 2013 (PEP8)
"""

import sys
import pylab

# if filename given as command line use this to read data,
# otherwise open 'data.txt'
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = "data.txt"

# Read data from text file (all rows must have same number of columns)
data = pylab.loadtxt(filename)

# Diagnostic output
print("Data in '%s' has %d rows and %d columns" %
      (filename, data.shape[0], data.shape[1]))

# plot columns 2 to : against column 1
x = data[:, 0]
ys = data[:, 1:]
pylab.plot(x, ys, '-o')
pylab.grid()
# if desired, save to file by uncommenting the next line(s)
# pylab.savefig('dataplot.png')
# pylab.savefig('dataplot.pdf')
pylab.show()

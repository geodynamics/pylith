#!/usr/bin/env python
import math
import numpy
import pdb
# pdb.set_trace()

# Get command-line arguments.
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--in_file", action="store", type="string",
                  dest="in_file", help="input file from GMT")
parser.add_option("-o", "--out_file", action="store", type="string",
                  dest="out_file", help="output VTK file")
parser.add_option("-e", "--elevation", action="store", type="float",
                  default=0.0, dest="elevation",
                  help="z-value assigned to segments [default: %default]")

(options, args) = parser.parse_args()

if (not options.in_file):
  parser.error("input file must be specified")

if (not options.out_file):
  parser.error("output file must be specified")

in_file = options.in_file
out_file = options.out_file
elevation = options.elevation

f = open(in_file, 'r')
g = open(out_file, 'w')

g.write('# vtk DataFile Version 2.0\n')
g.write('New Zealand coastlines.\n')
g.write('ASCII\n')
g.write('DATASET POLYDATA\n')

numvertices = 0
numlines = 0
vertsperline = []
totallineentries = 0
x = []
y = []
z = elevation

for line in f:
  if( line.startswith('#')):
    numlines += 1
    vertsperline.append(0)
    totallineentries += 1
  else:
    numvertices += 1
    vertsperline[numlines - 1] += 1
    totallineentries += 1
    data = line.split()
    x.append(float(data[0]))
    y.append(float(data[1]))

f.close()

g.write('POINTS %8d double\n' % numvertices)

for vertex in range(numvertices):
  g.write('%.10f' % x[vertex])
  g.write('  %.10f' % y[vertex])
  g.write('  %.10f' % z)
  g.write('\n')

g.write('LINES %8d  ' % numlines)
g.write('%8d\n' % totallineentries)

startvertex = 0
for line in range(numlines):
  g.write('%8d' % vertsperline[line])
  for vertex in range(startvertex, startvertex + vertsperline[line]):
    g.write('  %8d' % vertex)
  startvertex += vertsperline[line]

  g.write('\n')

g.close()
print "Number of vertices:  ", numvertices
print "Number of lines:  ", numlines
print "Number of line entries:  ", totallineentries

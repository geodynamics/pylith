#!/usr/bin/env python
"""
This script generates a set of line plots comparing the predicted solution to
the true solution.
"""

# The code requires the numpy, h5py, and matplotlib packages.
import numpy
import h5py
import matplotlib.pyplot as pyplot
# import pdb

# Define line colors and other parameters.
solnColor = "black"
predColors = ["cyan", "red", "blue", "green", "orange", "purple", "yellow"]
style = "lightbg"
fontsize = 8


def getSolution():
  """
  Function to get applied slip and fault coordinates.
  """

  # Open solution file and get slip and coordinates.
  solution = h5py.File(solutionFile, "r")
  solC = solution['geometry/vertices'][:]
  solV = solution['vertex_fields/slip'][:,:,0].flatten()

  # Sort by y-coordinate.
  solInds = numpy.argsort(solC[:,1])
  solCSort = solC[solInds,:]
  solVSort = solV[solInds]
  # Reverse sign, since this is right-lateral slip.
  solVSort *= -1.0
  solution.close()

  return (solCSort, solVSort)


# pdb.set_trace()

# The main part of the code is below.
# Get command-line arguments.
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--solution_file", action="store", type="string",
                  dest="solution_file",
                  help="HDF5 file with applied fault slip")
parser.add_option("-r", "--predicted_file", action="store", type="string",
                  dest="predicted_file", help="file with predicted solutions")

(options, args) = parser.parse_args()

if (not options.solution_file):
  parser.error("solution input file must be specified")

if (not options.predicted_file):
  parser.error("predicted input file must be specified")

solutionFile = options.solution_file
predictedFile = options.predicted_file

# Get applied fault slip.
(solCoords, solVals) = getSolution()

# Get predicted fault slip.
predicted = numpy.loadtxt(predictedFile, skiprows=1, dtype=numpy.float64)
predCoords = predicted[:,0:2]
predSlip = predicted[:, 2:]

# Generate figure, starting with true solution.
pyplot.plot(solCoords[:,1], solVals, color=solnColor)
pyplot.xlabel("Y-distance (m)")
pyplot.ylabel("Reverse slip (m)")

# Loop over predicted results, using a different line color for each.
numInv = predSlip.shape[1]

for inversion in range(numInv):
  pyplot.plot(predCoords[:,1], predSlip[:,inversion],
              color=predColors[inversion])


pyplot.show()

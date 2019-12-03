#!/usr/bin/env python
"""
This script generates a plot comparing the predicted solutions to
the true solution.
"""

# The code requires the numpy, h5py, and matplotlib packages.
import numpy
import h5py
import matplotlib.pyplot as pyplot

# Define line colors and other parameters.
solnColor = "black"
predColors = ["cyan", "red", "blue", "green", "orange", "purple", "yellow"]


# ----------------------------------------------------------------------
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

# ======================================================================
# The main part of the code is below.
# Get command-line arguments.
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", "--solution", action="store", type="string",
                  dest="solution_file",
                  help="HDF5 file with applied fault slip")
parser.add_option("-p", "--predicted", action="store", type="string",
                  dest="predicted_file", help="file with predicted solutions")

(options, args) = parser.parse_args()

if not options.solution_file:
  parser.error("Solution input file must be specified.")

if not options.predicted_file:
  parser.error("Predicted input file must be specified.")

solutionFile = options.solution_file
predictedFile = options.predicted_file

# Get applied fault slip.
(solCoords, solVals) = getSolution()

# Get predicted fault slip.
predicted = numpy.loadtxt(predictedFile)
predCoords = predicted[:,0:2]
predSlip = predicted[:, 2:]

# Generate figure, starting with true solution.
pyplot.plot(solCoords[:,1]/1.0e+3, solVals, linewidth=2, color=solnColor)
pyplot.xlabel("Distance Along Strike (km)")
pyplot.ylabel("Right-Lateral Slip (m)")

# Loop over predicted results, using a different line color for each.
numInv = predSlip.shape[1]

for inversion in range(numInv):
  pyplot.plot(predCoords[:,1]/1.0e+3, predSlip[:,inversion],
              color=predColors[inversion])


pyplot.show()
# pyplot.savefig("strikeslip_inversion.pdf")

#!/usr/bin/env nemesis
"""
This is an extremely simple example showing how to set up an inversion
using PyLith-generated Green's functions. In this simple example,
there are no data uncertainties, and for a penalty function we use a
simple minimum moment.
"""

# The code requires the numpy and h5py packages.
import numpy
import h5py


# ----------------------------------------------------------------------
def getImpResp():
  """
  Function to get impulse and response coordinates and values. Both are sorted
  along the y-dimension of the fault.
  """

  # Open impulse file and determine which fault vertices were used.
  impulses = h5py.File(impulseFile, "r", driver="sec2")
  impC = impulses['geometry/vertices'][:]
  impV = impulses['vertex_fields/slip'][:,:,0]
  impInds = numpy.nonzero(impV != 0.0)
  impCUsed = impC[impInds[1]]
  impVUsed = impV[impInds[0], impInds[1]]

  # Sort by y-coordinate.
  impInds = numpy.argsort(impCUsed[:,1])
  impCSort = impCUsed[impInds,:]
  impVSort = impVUsed[impInds]

  # Close impulse file and get responses.
  impulses.close()
  responses = h5py.File(responseFile, "r", driver="sec2")
  respC = responses['geometry/vertices'][:]
  respV = responses['vertex_fields/displacement'][:]
  respVIsort = respV[impInds,:,:]
  responses.close()

  return (impCSort, impVSort, respC, respVIsort)


# ----------------------------------------------------------------------
def getData():
  """
  Function to get data values and coordinates.
  """

  # Open data file. Since the response and data coordinates are in the same
  # order, we don't have to worry about sorting.
  data = h5py.File(dataFile, "r", driver="sec2")
  dataC = data['geometry/vertices'][:]
  dataV = data['vertex_fields/displacement'][:]
  data.close()

  return (dataC, dataV)


# ======================================================================
# The main part of the code is below.
# Get command-line arguments.
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--impulses", action="store", type="string",
                  dest="impulse_file", help="HDF5 file with fault GF info")
parser.add_option("-r", "--responses", action="store", type="string",
                  dest="response_file", help="HDF5 file with GF responses")
parser.add_option("-d", "--data", action="store", type="string",
                  dest="data_file", help="HDF5 file with data")
parser.add_option("-p", "--penalty", action="store", type="string",
                  dest="penalty_file", help="text file with penalty parameters")
parser.add_option("-o", "--output", action="store", type="string",
                  dest="output_file", help="text file with estimated slip")

(options, args) = parser.parse_args()

if not options.impulse_file:
  parser.error("Impulse input file must be specified.")

if not options.response_file:
  parser.error("Response input file must be specified.")

if not options.data_file:
  parser.error("Data input file must be specified.")

if not options.penalty_file:
  parser.error("penalty input file must be specified.")

if not options.output_file:
  parser.error("Output file must be specified.")

impulseFile = options.impulse_file
responseFile = options.response_file
dataFile = options.data_file
penaltyFile = options.penalty_file
outputFile = options.output_file

# Get GF info.
(impCoords, impVals, respCoords, respVals) = getImpResp()

# Get observed displacements and observation locations.
(dataCoords, dataVals) = getData()

# Get penalty parameters.
penalties = numpy.loadtxt(penaltyFile, dtype=numpy.float64)

# Determine matrix sizes and set up A-matrix.
numParams = impVals.shape[0]
numObs = 2 * dataVals.shape[1]
aMat = respVals.reshape((numParams, numObs)).transpose()

# Create diagonal matrix to use as the penalty.
parDiag = numpy.eye(numParams, dtype=numpy.float64)

# Data vector is a flattened version of the dataVals, plus the a priori
# values of the parameters (assumed to be zero).
dataVec = numpy.concatenate((dataVals.flatten(),
                             numpy.zeros(numParams, dtype=numpy.float64)))

# Determine number of inversions and create empty array to hold results.
numInv = penalties.shape[0]
invResults = numpy.zeros((numParams, 2 + numInv))
invResults[:,0:2] = impCoords
head = "# X_Coord Y_Coord"

# Loop over number of inversions.
for inversion in range(numInv):
  penalty = penalties[inversion]
  head += " Penalty=%g" % penalty

  # Scale diagonal by penalty parameter, and stack A-matrix with penalty matrix.
  penMat = penalty * parDiag
  designMat = numpy.vstack((aMat, penMat))
  designMatTrans = designMat.transpose()

  # Form generalized inverse matrix.
  genInv = numpy.dot(numpy.linalg.inv(numpy.dot(designMatTrans, designMat)),
                                      designMatTrans)

  # Solution is matrix-vector product of generalized inverse with data vector.
  solution = numpy.dot(genInv, dataVec)
  invResults[:,2 + inversion] = solution

  # Compute predicted results and residual.
  predicted = numpy.dot(aMat, solution)
  residual = dataVals.flatten() - predicted
  residualNorm = numpy.linalg.norm(residual)
  print "Penalty parameter:  %g" % penalty
  print "  Residual norm:    %g" % residualNorm

head += "\n"

# Output results.
f = open(outputFile, "w")
f.write(head)
numpy.savetxt(f, invResults, fmt="%14.6e")
f.close()

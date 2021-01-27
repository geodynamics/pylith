#!/usr/bin/env python
# -*- Python -*- (syntax highlighting)
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

# This script creates a log-log plot of weighted data misfit vs. penalty misfit.
# The code requires the numpy and matplotlib packages.

import numpy
import matplotlib.pyplot as pyplot
import math

# Define line colors and other parameters.
lineColor = "blue"

# ----------------------------------------------------------------------
def readInversionSummary():
  """
  Function to read inversion results from text file.
  """

  # Open inversion summary file and get misfits.
  data = numpy.loadtxt(summaryFile, dtype=numpy.float64)
  dataWeightResid = data[:, 2]
  penaltyResid = data[:, 3]

  # Sort by penalty residual.
  inds = numpy.argsort(penaltyResid)
  dataWeightResidSort = dataWeightResid[inds]
  penaltyResidSort = penaltyResid[inds]

  return (dataWeightResidSort, penaltyResidSort)


# ======================================================================
# The main part of the code is below.
# Get command-line arguments.
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-s", "--summary", action="store", type="string",
                  dest="summary_file",
                  help="Text file with inversion summary")

(options, args) = parser.parse_args()

if not options.summary_file:
  parser.error("Summary input file must be specified.")

summaryFile = options.summary_file

# Get misfits.
(dataWeightResid, penaltyResid) = readInversionSummary()

# Generate figure, starting with true solution.
pyplot.loglog(penaltyResid, dataWeightResid, linewidth=2, color=lineColor,
              marker="o")
pyplot.xlabel("Log(Penalty Residual)")
pyplot.ylabel("Log(Weighted Data Residual)")
pyplot.xlim(left=10.0, right=30.0)

pyplot.show()
# pyplot.savefig("reverse_inversion.pdf")

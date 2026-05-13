#!/usr/bin/env nemesis
# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information.
# =================================================================================================
"""This script creates a log-log plot of weighted data misfit vs. penalty misfit.

The code requires the numpy and matplotlib packages.
"""

import numpy
import matplotlib.pyplot as pyplot

SUMMARY_FILENAME = "output/step07-inversion-summary.txt"


# ----------------------------------------------------------------------
def readInversionSummary():
    """Function to read inversion results from text file."""

    # Open inversion summary file and get misfits.
    data = numpy.loadtxt(SUMMARY_FILENAME)
    dataWeightResidual = data[:, 2]
    penaltyResidual = data[:, 3]

    # Sort by penalty residual.
    inds = numpy.argsort(penaltyResidual)
    dataWeightResidSort = dataWeightResidual[inds]
    penaltyResidSort = penaltyResidual[inds]

    return (dataWeightResidSort, penaltyResidSort)


# Get misfits.
dataWeightResidual, penaltyResidual = readInversionSummary()

# Generate figure, starting with true solution.
pyplot.loglog(penaltyResidual, dataWeightResidual, linewidth=2, marker="o")
pyplot.xlabel("log(Penalty Residual)")
pyplot.ylabel("log(Weighted Data Residual)")
# pyplot.xlim(left=10.0, right=30.0)

pyplot.show()
# pyplot.savefig("reverse_inversion.pdf")

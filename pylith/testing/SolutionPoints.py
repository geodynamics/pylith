# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

# @file pylith/testing/SolutionPoints.py
##
# @brief Check displacement solution point output from PyLith.

import numpy
import h5py


def check_displacements(testcase, filename, npoints, spaceDim):
    """Check displacements.
    """
    h5 = h5py.File(filename, "r", driver="sec2")

    # Check vertices
    vertices = h5['geometry/vertices'][:]
    (nvertices, spaceDim) = vertices.shape
    testcase.assertEqual(npoints, nvertices)
    testcase.assertEqual(spaceDim, spaceDim)

    # Check displacement solution
    toleranceAbsMask = 0.1
    tolerance = 1.0e-5

    dispE = testcase.calcDisplacementPoints(vertices)
    disp = h5['vertex_fields/displacement'][:]

    (nstepsE, nverticesE, ncompsE) = dispE.shape
    (nsteps, nvertices, ncomps) = disp.shape
    testcase.assertEqual(nstepsE, nsteps)
    testcase.assertEqual(nverticesE, nvertices)
    testcase.assertEqual(ncompsE, ncomps)

    for istep in range(nsteps):
        for icomp in range(ncomps):
            okay = numpy.zeros((nvertices,), dtype=numpy.bool)

            maskR = numpy.abs(dispE[istep, :, icomp]) > toleranceAbsMask
            ratio = numpy.abs(
                1.0 - disp[istep, maskR, icomp] / dispE[istep, maskR, icomp])
            if len(ratio) > 0:
                okay[maskR] = ratio < tolerance

            maskD = ~maskR
            diff = numpy.abs(disp[istep, maskD, icomp] -
                             dispE[istep, maskD, icomp])
            if len(diff) > 0:
                okay[maskD] = diff < tolerance

            if numpy.sum(okay) != nvertices:
                print("Error in component %d of displacement field at time step %d." % (
                    icomp, istep))
                print("Expected values: ", dispE[istep, :, :])
                print("Output values: ", disp[istep, :, :])
                print("Expected values (not okay): ",
                      dispE[istep, ~okay, icomp])
                print("Computed values (not okay): ",
                      disp[istep, ~okay, icomp])
                print("Relative diff (not okay): ", diff[~okay])
                print("Coordinates (not okay): ", vertices[~okay, :])
                h5.close()
            testcase.assertEqual(nvertices, numpy.sum(okay))

    h5.close()
    return


def check_stations(testcase, filename, stationsE):
    """Check station names.
    """
    h5 = h5py.File(filename, "r", driver="sec2")
    stations = h5['stations'][:]
    h5.close()

    nstationsE = len(stations)
    nstations = stations.shape[0]
    testcase.assertEqual(nstationsE, nstations)

    for stE, st in zip(stationsE, stations):
        testcase.assertEqual(stE, st)

    return

# End of file

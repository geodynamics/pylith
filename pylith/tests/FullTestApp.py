#!/usr/bin/env nemesis
#
# ======================================================================
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
# ======================================================================
#

# @file pylith/tests/FullTestApp.py

# @brief Python application for Python full-scale tests.

import unittest
import numpy

from pylith.tests import has_h5py
from pylith.apps.PyLithApp import PyLithApp


# ----------------------------------------------------------------------------------------------------------------------
class TestCase(unittest.TestCase):
    """
    Generic test case for full-scale test.
    """
    NAME = None  # Set in child class.
    VERBOSITY = 0
    RUN_PYLITH = True

    def setUp(self):
        """
        Setup for test.
        """
        return

    def run_pylith(self, testName, args, generatedb=None):
        if self.RUN_PYLITH:
            if self.VERBOSITY > 0:
                print("Running Pylith with args '{}' ...".format(" ".join(args)))
            run_pylith(testName, args, generatedb)
        return

    @staticmethod
    def parse_args():
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--verbose", action="store_true", dest="verbosity", default=0)
        parser.add_argument("--skip-pylith-run", action="store_false", dest="run_pylith", default=True)
        args = parser.parse_args()
        TestCase.VERBOSITY = args.verbosity
        TestCase.RUN_PYLITH = args.run_pylith
        return


# ----------------------------------------------------------------------------------------------------------------------
class Example(TestCase):
    """
    Base class for running an example.

    Need one test_* method to instantiate object and run PyLith via setUp().
    """
    NAME = None
    PYLITH_ARGS = None

    def setUp(self):
        TestCase.setUp(self)
        TestCase.run_pylith(self, self.NAME, self.PYLITH_ARGS)

    def test_example(self):
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestDriver(object):
    """
    Driver application for running full-scale tests.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self):
        """
        Constructor.
        """
        return

    def main(self):
        """
        Run the test suite.
        """
        success = unittest.TextTestRunner(verbosity=2).run(self._suite()).wasSuccessful()

        if not success:
            import sys
            sys.exit(1)
        return


# ----------------------------------------------------------------------------------------------------------------------
class HDF5Checker(object):

    def __init__(self, filename, testcase, mesh, ratio_tolerance=1.e-5, diff_tolerance=1e-5):
        """Constructor.
        """
        import h5py
        self.h5 = h5py.File(filename, "r")
        self.testcase = testcase
        self.exactsoln = testcase.exactsoln
        self.mesh = mesh
        self.ratio_tolerance = ratio_tolerance
        self.diff_tolerance = diff_tolerance
        self.vertices = None
        self.cellCentroids = None
        return

    def _getVertices(self):
        """Get vertices, reading from file if necessary.
        """
        if self.vertices is None:
            self.testcase.assertTrue("geometry" in self.h5.keys())
            self.testcase.assertTrue("vertices" in self.h5["geometry"].keys())
            vertices = self.h5["geometry/vertices"][:]
            (nvertices, spaceDim) = vertices.shape
            self.testcase.assertEqual(self.mesh['nvertices'], nvertices)
            self.testcase.assertEqual(self.exactsoln.SPACE_DIM, spaceDim)
            self.vertices = vertices
        return self.vertices

    def _getCellCentroids(self):
        """Get cell centroids, reading cells and vertices from file if necessary.
        """
        if self.cellCentroids is None:
            vertices = self._getVertices()
            self.testcase.assertTrue("topology" in self.h5.keys())
            self.testcase.assertTrue("cells" in self.h5["topology"].keys())
            cells = self.h5["topology/cells"][:].astype(numpy.int)
            ncells, ncorners = cells.shape
            centroids = numpy.zeros((ncells, self.exactsoln.SPACE_DIM), dtype=numpy.float64)
            for icorner in range(ncorners):
                centroids[:,:] += vertices[cells[:, icorner],:]
            centroids /= float(ncorners)
            self.cellCentroids = centroids
        return self.cellCentroids

    def checkVertexField(self, fieldName):
        self.testcase.assertTrue("vertex_fields" in self.h5.keys())
        fieldsH5 = self.h5["vertex_fields"].keys()
        self.testcase.assertTrue(fieldName in fieldsH5,
                                 "Could not find field '{}' in vertex fields {}".format(fieldName, fieldsH5))
        field = self.h5["vertex_fields/" + fieldName][:]

        vertices = self._getVertices()
        (nvertices, spaceDim) = vertices.shape
        fieldE = self.exactsoln.getField(fieldName, vertices)
        mask = None
        if "getMask" in dir(self.exactsoln):
            mask = self.exactsoln.getMask(fieldName, vertices)
        self._checkField(fieldName, fieldE, field, vertices, mask)
        return

    def checkCellField(self, fieldName):
        self.testcase.assertTrue("cell_fields" in self.h5.keys(),
                                 "Missing 'cell_fields'. Groups: {}".format(self.h5.keys()))
        fieldsH5 = self.h5["cell_fields"].keys()
        self.testcase.assertTrue(fieldName in fieldsH5,
                                 "Could not find field '{}' in cell fields {}".format(fieldName, fieldsH5))
        field = self.h5["cell_fields/" + fieldName][:]

        centroids = self._getCellCentroids()
        (ncells, spaceDim) = centroids.shape
        fieldE = self.exactsoln.getField(fieldName, centroids)
        mask = None
        self._checkField(fieldName, fieldE, field, centroids, mask)
        return

    def _checkField(self, fieldName, fieldE, field, pts, maskField):
        (nstepsE, nptsE, ncompsE) = fieldE.shape
        (nsteps, npts, ncomps) = field.shape
        self.testcase.assertEqual(nstepsE, nsteps, msg="Expected {} time steps, got {} for field {}".format(
            nstepsE, nsteps, fieldName))
        self.testcase.assertEqual(nptsE, npts, msg="Expected {} points, got {} for field {}".format(
            nptsE, npts, fieldName))
        self.testcase.assertEqual(ncompsE, ncomps, msg="Expected {} components, got {} for field: {}".format(
            ncompsE, ncomps, fieldName))

        toleranceAbsMask = 0.1
        ratio_tolerance = self.ratio_tolerance
        diff_tolerance = self.diff_tolerance
        #maskZero = fieldE != 0.0
        maskZero = numpy.abs(fieldE) > 1e-15
        scale = numpy.mean(numpy.abs(fieldE[maskZero].ravel())) if numpy.sum(maskZero) > 0 else 1.0
        for istep in range(nsteps):
            for icomp in range(ncomps):
                okay = numpy.zeros((npts,), dtype=numpy.bool)
                ratio = numpy.zeros((npts,))

                maskR = numpy.abs(fieldE[istep,:, icomp]) > toleranceAbsMask
                ratio[maskR] = numpy.abs(1.0 - field[istep, maskR, icomp] / fieldE[istep, maskR, icomp])
                if numpy.sum(maskR) > 0:
                    okay[maskR] = ratio[maskR] < ratio_tolerance

                maskD = ~maskR
                diff = numpy.abs(field[istep,:, icomp] - fieldE[istep,:, icomp]) / scale
                if numpy.sum(maskD) > 0:
                    okay[maskD] = diff[maskD] < diff_tolerance

                if not maskField is None:
                    okay[maskField] = True

                if numpy.sum(okay) != npts:
                    print("Error in component {} of field '{}' at time step {}.".format(icomp, fieldName, istep))
                    # Debug Output
                #    print("Expected values: ", fieldE[istep, :, :])
                #    print("Output values: ", field[istep, :, :])

                    print("Total # not okay: %d" % numpy.sum(~okay))
                    n_okay_maskR = numpy.logical_and(~okay, maskR)
                    if numpy.sum(n_okay_maskR) > 0:
                        print("Ratio (maskR), # not okay: %d" % numpy.sum(n_okay_maskR))
                        print("Expected values (not okay): %s" % fieldE[istep, n_okay_maskR, icomp])
                        print("Computed values (not okay): %s" % field[istep, n_okay_maskR, icomp])
                        print("Ratio (not okay): %s" % ratio[n_okay_maskR])
                        print("Ratio Coordinates (not okay): %s" % pts[n_okay_maskR, :])
                        print("Tolerance Absolute Mask: %s" % toleranceAbsMask)
                        print("Ratio Tolerance: %s" % ratio_tolerance)

                    n_okay_maskD = numpy.logical_and(~okay, maskD)
                    if numpy.sum(n_okay_maskD) > 0:
                        print("Relative Diff (maskD), # not okay: %d" % numpy.sum(n_okay_maskD))
                        print("Expected values (not okay): %s" % fieldE[istep, n_okay_maskD, icomp])
                        print("Computed values (not okay): %s" % field[istep, n_okay_maskD, icomp])
                        print("Relative diff (not okay): %s" % diff[n_okay_maskD])
                        print("Relative diff Coordinates (not okay): %s" % pts[n_okay_maskD, :])
                        print("Diff Tolerance: %f" % diff_tolerance)
                        print("Scale: %10.4e" % scale)
                    self.testcase.assertEqual(npts, numpy.sum(okay))

        return


# ----------------------------------------------------------------------------------------------------------------------
def check_data(filename, testcase, mesh, vertexFields=[], cellFields=[], ratio_tolerance=1.e-5, diff_tolerance=1.e-5):
    """Check vertex and cell fields in specified file.
    """
    separateChecker = False
    if not has_h5py():
        return

    if type(ratio_tolerance) is dict:
        separateChecker = True
        defaultRatio = 1e-5

    checker = HDF5Checker(filename, testcase, mesh, ratio_tolerance=ratio_tolerance, diff_tolerance=diff_tolerance)
    for field in vertexFields:
        if separateChecker:
            checker = HDF5Checker(filename, testcase, mesh, ratio_tolerance=ratio_tolerance.get(str(field), defaultRatio), diff_tolerance=diff_tolerance.get(str(field), defaultRatio))
        if testcase.VERBOSITY > 0:
            print("Checking vertex field '{}' in file {}.".format(field, filename))
        checker.checkVertexField(field)
    for field in cellFields:
        if separateChecker:
            checker = HDF5Checker(filename, testcase, mesh, ratio_tolerance=ratio_tolerance.get(str(field), defaultRatio), diff_tolerance=diff_tolerance.get(str(field), defaultRatio))
        if testcase.VERBOSITY > 0:
            print("Checking cell field '{}' in file {}.".format(field, filename))
        checker.checkCellField(field)
    return


# ----------------------------------------------------------------------------------------------------------------------
def run_pylith(appName, cfgfiles=[], dbClass=None, nprocs=1):
    """
    Helper function to generate spatial databases and run PyLith.
    """
    # Skip running if already run.
    if str(appName) in dir(run_pylith):
        return

    # Generate spatial databases if necessary.
    if not dbClass is None:
        db = dbClass()
        db.run()

    # Limit number of processes to number of local CPUs or maximum specified by environment.
    import os
    if "MAX_PYLITH_PROCS" in os.environ:
        appNumProcs = min(int(os.environ["MAX_PYLITH_PROCS"]), nprocs)
        if appNumProcs < nprocs:
            print("WARNING: Detected environment with MAX_PYLITH_PROCS=%d. Reducing number of processes from %d to %d." % (
                appNumProcs, nprocs, appNumProcs))
    else:
        import pylith.utils.CollectVersionInfo
        import multiprocessing

        cpuCount = multiprocessing.cpu_count()
        mpiVersion = pylith.utils.CollectVersionInfo.CollectVersionInfo._collectVersionMPI()
        if mpiVersion["implementation"] == "OpenMPI" and mpiVersion["standard"].startswith("3"):
            cpuCount /= 2  # Assume hyperthreading is turned on and OpenMPI 3 doesn't allow oversubscribing

        appNumProcs = min(cpuCount, nprocs)
        if appNumProcs < nprocs:
            print("WARNING: Detected %d CPUs. Reducing number of processes from %d to %d." %
                  (appNumProcs, nprocs, appNumProcs))

    # Run Pylith
    app = PyLithApp()
    app.nodes = appNumProcs
    setattr(run_pylith, str(appName), True)
    app.run(argv=["pylith"] + cfgfiles)
    return


# End of file

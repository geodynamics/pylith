#!/usr/bin/env nemesis
#
# ======================================================================
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
# ======================================================================
#

# @file pylith/testing/FullTestApp.py

# @brief Python application for Python full-scale tests.

import unittest
import numpy

from pylith.testing import has_h5py
from pylith.apps.PyLithApp import PyLithApp


# -------------------------------------------------------------------------------------------------
class MeshEntity(object):

    def __init__(self, ncells, ncorners, nvertices):
        """Mesh entity (domain, boundary, material, fault, etc).

        Args:
            ncells (int)
                Number of cells in entity.
            ncorners (int)
                Number of vertices in cell.
            nvertices (int)
                Number of vertices in entity.
        """
        self.ncells = ncells
        self.ncorners = ncorners
        self.nvertices = nvertices

# -------------------------------------------------------------------------------------------------
class Check(object):

    def __init__(self, mesh_entities, filename=None, mesh=None, vertex_fields=[], cell_fields=[], exact_soln=None, tolerance=None, defaults={}):
        """Set parameters for checking PyLith output.

        Args
            mesh_entities (list of str):
                List of mesh entities corresponding to names of domain, materials, boundaries, etc.
            filename (str):
                File name template with name and mesh_entity keys.
            mesh (object):
                Object with number of points for mesh entities.
            vertex_fields (list of str):
                List of vertex fields to check.
            cell_fields (list of str):
                List of cell fields to check.
            exact_soln (object):
                Object with functions returning expected values for fields.
            defaults (dict):
                Dictionary with default values.
        """
        self.tolerance = 1.0e-5
        self.zero_tolerance = 1.0e-10
        self.vertex_fields = []
        self.cell_fields = []

        for key, value in defaults.items():
            setattr(self, key, value)

        self.mesh_entities = mesh_entities

        if filename:
            self.filename = filename
        if mesh:
            self.mesh = mesh
        if vertex_fields:
            self.vertex_fields = vertex_fields
        if cell_fields:
            self.cell_fields = cell_fields
        if exact_soln:
            self.exact_soln = exact_soln
        if tolerance:
            self.tolerance = tolerance

# -------------------------------------------------------------------------------------------------
class FullTestCase(unittest.TestCase):
    """Generic test case for full-scale test.
    """
    VERBOSITY = 0
    RUN_PYLITH = True

    def setUp(self):
        super().setUp()
        self.name = None
        self.mesh = None

    def test_output(self):
        for check in self.checks:
            for mesh_entity in check.mesh_entities:
                filename = check.filename.format(name=self.name, mesh_entity=mesh_entity)
                with self.subTest(filename=filename):
                    check_data(self, filename, check, mesh_entity, check.mesh.ENTITIES[mesh_entity])

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
        FullTestCase.VERBOSITY = args.verbosity
        FullTestCase.RUN_PYLITH = args.run_pylith
        return


# ----------------------------------------------------------------------------------------------------------------------
class TestDriver(object):
    """Driver application for running full-scale tests.
    """

    # PUBLIC METHODS /////////////////////////////////////////////////////

    def __init__(self):
        """Constructor.
        """
        return

    def main(self):
        """Run the test suite.
        """
        success = unittest.TextTestRunner(verbosity=2).run(self._suite()).wasSuccessful()

        if not success:
            import sys
            sys.exit(1)
        return


# ----------------------------------------------------------------------------------------------------------------------
class HDF5Checker(object):

    def __init__(self, filename, testcase, mesh):
        """Constructor.
        """
        import h5py
        self.debug = False

        self.h5 = h5py.File(filename, "r")
        self.testcase = testcase
        self.mesh = mesh
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
            self.testcase.assertEqual(self.mesh.nvertices, nvertices)
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
            (nvertices, spaceDim) = vertices.shape
            centroids = numpy.zeros((ncells, spaceDim), dtype=numpy.float64)
            for icorner in range(ncorners):
                centroids[:,:] += vertices[cells[:, icorner],:]
            centroids /= float(ncorners)
            self.cellCentroids = centroids
        return self.cellCentroids

    def checkVertexField(self, fieldName, mesh_entity, exact_soln, tolerance, zero_tolerance):
        self.testcase.assertTrue("vertex_fields" in self.h5.keys())
        fieldsH5 = self.h5["vertex_fields"].keys()
        self.testcase.assertTrue(fieldName in fieldsH5,
                                 f"Could not find field in vertex fields {fieldsH5}")
        field = self.h5["vertex_fields/" + fieldName][:]

        vertices = self._getVertices()
        (nvertices, spaceDim) = vertices.shape
        self.testcase.assertEqual(exact_soln.SPACE_DIM, spaceDim)
        fieldE = exact_soln.getField(fieldName, mesh_entity, vertices)
        mask = None
        if "getMask" in dir(exact_soln):
            mask = exact_soln.getMask(fieldName, mesh_entity, vertices)
        self._checkField(fieldE, field, mask, tolerance, zero_tolerance)
        return

    def checkCellField(self, fieldName, mesh_entity, exact_soln, tolerance, zero_tolerance):
        self.testcase.assertTrue("cell_fields" in self.h5.keys(),
                                 f"Missing 'cell_fields'. Groups: {self.h5.keys()}".format)
        fieldsH5 = self.h5["cell_fields"].keys()
        self.testcase.assertTrue(fieldName in fieldsH5,
                                 f"Could not find field in cell fields {fieldsH5}")
        field = self.h5["cell_fields/" + fieldName][:]

        centroids = self._getCellCentroids()
        fieldE = exact_soln.getField(fieldName, mesh_entity, centroids)
        mask = None
        self._checkField(fieldE, field, mask, tolerance, zero_tolerance)
        return

    def _checkField(self, fieldE, field, maskField, tolerance, zero_tolerance):
        (nstepsE, nptsE, ncompsE) = fieldE.shape
        (nsteps, npts, ncomps) = field.shape
        self.testcase.assertEqual(nstepsE, nsteps, msg="Mismatch in number of time steps")
        self.testcase.assertEqual(nptsE, npts, msg="Mismatch in number of points")
        self.testcase.assertEqual(ncompsE, ncomps, msg="Mismatch in number of components")

        scale = numpy.mean(numpy.abs(fieldE).ravel())
        vtolerance = scale*tolerance if scale > zero_tolerance else tolerance
        okay = numpy.abs(field - fieldE) < vtolerance

        if not maskField is None:
            okay[:, maskField, :] = True

        msg = []
        if numpy.sum(okay) != nsteps*npts*ncomps:
            if self.debug:
                msg += [f"Expected values: {fieldE[:]}"]
                msg += [f"Computed values: {field[:]}"]

            msg += [""]
            msg += [f"Expected values (not okay): {fieldE[~okay]}"]
            msg += [f"Computed values (not okay): {field[~okay]}"]
            #msg += [f"Coordinates: {pts[~okay, :]}"]
            msg += [f"Tolerance: {vtolerance}"]
        self.testcase.assertEqual(nsteps*npts*ncomps, numpy.sum(okay), msg="\n".join(msg))


# ----------------------------------------------------------------------------------------------------------------------
def check_data(testcase, filename, check, mesh_entity, mesh):
    """Check vertex and cell fields in specified file.
    """
    if not has_h5py():
        return

    checker = HDF5Checker(filename, testcase, mesh)

    for field in check.vertex_fields:
        field_tolerance = check.tolerance[field] if isinstance(check.tolerance, dict) else check.tolerance
        field_zero_tolerance = check.zero_tolerance[field] if isinstance(check.zero_tolerance, dict) else check.zero_tolerance
        with testcase.subTest(vertex_field=field):
            checker.checkVertexField(field, mesh_entity, check.exact_soln, field_tolerance, field_zero_tolerance)

    for field in check.cell_fields:
        field_tolerance = check.tolerance[field] if isinstance(check.tolerance, dict) else check.tolerance
        field_zero_tolerance = check.zero_tolerance[field] if isinstance(check.zero_tolerance, dict) else check.zero_tolerance
        with testcase.subTest(cell_field=field):
            checker.checkCellField(field, mesh_entity, check.exact_soln, field_tolerance, field_zero_tolerance)
    return


# ----------------------------------------------------------------------------------------------------------------------
def run_pylith(appName, cfgfiles=[], dbClass=None, nprocs=1):
    """Helper function to generate spatial databases and run PyLith.
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

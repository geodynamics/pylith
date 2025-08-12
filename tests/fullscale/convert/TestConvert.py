# =================================================================================================
# This code is part of PyLith, developed through the Computational Infrastructure
# for Geodynamics (https://github.com/geodynamics/pylith).
#
# Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
# All rights reserved.
#
# See https://mit-license.org/ and LICENSE.md and for license information. 
# =================================================================================================

import unittest
import numpy
import h5py

from pylith.apps.ConvertMeshApp import ConvertMeshApp


class TestConvert(unittest.TestCase):
    """Test converting mesh format.
    """
    def setUp(self):
        """Setup for test.
        """
        run_convert("quad", ["convert_quad.cfg"])

    def test_quad(self):
        """Check mesh.
        """
        self._checkMeshHDF5("quad_small.h5", dim=2, numCells=12, numCorners=4, numVertices=17)

    def _checkMeshHDF5(self, filename, dim, numCells, numCorners, numVertices):
        groups = (
            "/topologies",
            "/topologies/domain",
            "/topologies/domain/distributions",
            "/topologies/domain/dms",
            "/topologies/domain/dms/coordinateDM",
            "/topologies/domain/dms/coordinateDM/section",
            "/topologies/domain/dms/coordinateDM/vecs",
            "/topologies/domain/dms/coordinateDM/vecs/coordinates",
            "/topologies/domain/labels",
            "/topologies/domain/labels/celltype",
            "/topologies/domain/labels/fault",
            "/topologies/domain/labels/material-id",
            "/topologies/domain/topology",
            "/topologies/domain/topology/strata",
            "/viz",
            "/viz/topology",
        )
        h5 = h5py.File(filename)
        for group in groups:
            item = h5.get(group)
            self.assertTrue(isinstance(item, h5py.Group))

        cells = h5["/viz/topology/cells"]
        self.assertEqual(dim, cells.attrs["cell_dim"])
        self.assertEqual(numCorners, cells.attrs["cell_corners"])
        self.assertEqual(numCells, cells[:].shape[0])
        self.assertEqual(numCorners, cells[:].shape[1])

        coords = h5["/topologies/domain/dms/coordinateDM/vecs/coordinates/coordinates"][:]
        self.assertEqual(numVertices*dim, coords.shape[0])

        cone_sizes = h5["/topologies/domain/topology/strata/2/cone_sizes"][:]
        self.assertEqual(numCells, cone_sizes[0, 0])
        self.assertEqual(numCorners, cone_sizes[0, 2])


# ----------------------------------------------------------------------
def run_convert(appName, cfgfiles):
    """Helper function to run pylith_convert.
    """
    if str(appName) in dir(run_convert):
        return

    app = ConvertMeshApp()
    setattr(run_convert, str(appName), True)
    app.run(argv=["test_convert"] + cfgfiles)


# End of file

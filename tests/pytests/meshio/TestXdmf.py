#!/usr/bin/env python
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

## @file tests/pytests/meshio/TestXdmf.py

## @brief Unit testing of Python Xdmf object.

import unittest

from pylith.meshio.Xdmf import Xdmf

# ----------------------------------------------------------------------
class TestXdmf(unittest.TestCase):
    """Unit testing of Python Xdmf object.
    """

    def test_constructor(self):
        """Test constructor.
        """
        xdmf = Xdmf()
        return

      
    def test_write(self):
        FILES = [
          "tri3.h5",
          "tri3_vertex.h5",
          "tri3_cell.h5",
          "tri3_points.h5",
          "tri3_points_vertex.h5",
          "tri3_surf.h5",
          "tri3_surf_vertex.h5",
          "tri3_surf_cell.h5",
          "quad4.h5",
          "quad4_vertex.h5",
          "quad4_cell.h5",
          "quad4_points.h5",
          "quad4_points_vertex.h5",
          "quad4_surf.h5",
          "quad4_surf_vertex.h5",
          "quad4_surf_cell.h5",
          "tet4.h5",
          "tet4_vertex.h5",
          "tet4_cell.h5",
          "tet4_points.h5",
          "tet4_points_vertex.h5",
          "tet4_surf.h5",
          "tet4_surf_vertex.h5",
          "tet4_surf_cell.h5",
          "hex8.h5",
          "hex8_vertex.h5",
          "hex8_cell.h5",
          "hex8_points.h5",
          "hex8_points_vertex.h5",
          "hex8_surf.h5",
          "hex8_surf_vertex.h5",
          "hex8_surf_cell.h5",
        ]

        import os

        xdmf = Xdmf()
        prefixData = "data"
        prefix = ""
        if not os.path.isdir("data"):
            prefixData = os.path.join("meshio", "data")
            prefix = "meshio"
        for filenameH5 in FILES:
            pathH5 = os.path.join(prefixData, filenameH5)
            filenameXdmf = os.path.join(prefix, filenameH5.replace(".h5", ".xmf"))
            xdmf.write(pathH5, filenameXdmf, verbose=False)

            filenameXdmfE = pathH5.replace(".h5", ".xmf")
            self._check(filenameXdmfE, filenameXdmf)


    def _check(self, filenameE, filename):
        fin = open(filenameE, "r")
        linesE = fin.readlines()
        fin.close()
        
        fin = open(filename, "r")
        lines = fin.readlines()
        fin.close()

        self.assertEqual(len(linesE), len(lines), "Number of lines for Xdmf file '%s' doesn't match." % filename)

        numLines = len(linesE)
        for i in range(numLines):
            self.assertEqual(linesE[i], lines[i], "Line %d of file '%s' doesn't match." % (i, filename))
        return
  

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestXdmf))
    unittest.TextTestRunner(verbosity=2).run(suite)


# End of file 

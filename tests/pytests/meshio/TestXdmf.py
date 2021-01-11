#!/usr/bin/env python
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

## @file tests/pytests/meshio/TestXdmf.py

## @brief Unit testing of Python Xdmf object.

import unittest

from pylith.meshio.Xdmf import Xdmf

# ----------------------------------------------------------------------
class TestXdmf(unittest.TestCase):
    """
    Unit testing of Python Xdmf object.
    """

    def test_constructor(self):
        """
       Test constructor.
        """
        xdmf = Xdmf()
        return

      
    def test_write(self):
        files = [
          "data/tri3.h5",
          "data/tri3_vertex.h5",
          "data/tri3_cell.h5",
          "data/tri3_points.h5",
          "data/tri3_points_vertex.h5",
          "data/tri3_surf.h5",
          "data/tri3_surf_vertex.h5",
          "data/tri3_surf_cell.h5",
          "data/quad4.h5",
          "data/quad4_vertex.h5",
          "data/quad4_cell.h5",
          "data/quad4_points.h5",
          "data/quad4_points_vertex.h5",
          "data/quad4_surf.h5",
          "data/quad4_surf_vertex.h5",
          "data/quad4_surf_cell.h5",
          "data/tet4.h5",
          "data/tet4_vertex.h5",
          "data/tet4_cell.h5",
          "data/tet4_points.h5",
          "data/tet4_points_vertex.h5",
          "data/tet4_surf.h5",
          "data/tet4_surf_vertex.h5",
          "data/tet4_surf_cell.h5",
          "data/hex8.h5",
          "data/hex8_vertex.h5",
          "data/hex8_cell.h5",
          "data/hex8_points.h5",
          "data/hex8_points_vertex.h5",
          "data/hex8_surf.h5",
          "data/hex8_surf_vertex.h5",
          "data/hex8_surf_cell.h5",
        ]

        import os

        xdmf = Xdmf()
        for filenameH5 in files:
            filenameXdmf = os.path.split(filenameH5)[-1].replace(".h5", ".xmf")
            xdmf.write(filenameH5, filenameXdmf, verbose=False)

            filenameXdmfE = "data/" + filenameXdmf
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
  

# End of file 

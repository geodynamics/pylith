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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

## @file unittests/pytests/topology/TestMeshImporter.py

## @brief Unit testing of MeshImporter object.

import unittest

from pylith.topology.MeshImporter import MeshImporter

# ----------------------------------------------------------------------
class TestMeshImporter(unittest.TestCase):
  """
  Unit testing of MeshIO object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    importer = MeshImporter()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.topology.MeshImporter import mesh_generator
    g = mesh_generator()
    return


# End of file 

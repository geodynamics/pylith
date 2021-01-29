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

## @file tests/pytests/faults/TestTractPerturbation.py

## @brief Unit testing of TractPerturbation object.

import unittest

from pylith.faults.TractPerturbation import TractPerturbation

# ----------------------------------------------------------------------
class TestTractPerturbation(unittest.TestCase):
  """
  Unit testing of TractPerturbation object.
  """

  def test_constructor(self):
    """
    Test constructor.
    """
    tract = TractPerturbation()
    return


  def test_configure(self):
    """
    Test initialize().
    """
    from spatialdata.spatialdb.SimpleDB import SimpleDB
    from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
    from pythia.pyre.units.time import second

    ioInitial = SimpleIOAscii()
    ioInitial.inventory.filename = "tri3_initialtractions.spatialdb"
    ioInitial._configure()
    dbInitial = SimpleDB()
    dbInitial.inventory.iohandler = ioInitial
    dbInitial.inventory.label = "initial tractions"
    dbInitial._configure()
    
    ioChange = SimpleIOAscii()
    ioChange.inventory.filename = "tri3_changetractions.spatialdb"
    ioChange._configure()
    dbChange = SimpleDB()
    dbChange.inventory.iohandler = ioChange
    dbChange.inventory.label = "traction change"
    dbChange._configure()
    
    tract = TractPerturbation()
    tract.inventory.dbInitial = dbInitial
    tract.inventory.dbChange = dbChange
    tract._configure()
    return


  def test_factory(self):
    """
    Test factory method.
    """
    from pylith.faults.TractPerturbation import traction_perturbation
    fn = traction_perturbation()
    return


# End of file 

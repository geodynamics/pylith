#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/problems/SubfieldPressure.py
##
## @brief Python object for pressure subfield.
##
## Factory: subfield.

from pylith.topology.Subfield import Subfield

# SubfieldPressure class
class SubfieldPressure(Subfield):
  """
  Python object for pressure subfield.

  Factory: subfield.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="subfieldpressure"):
    """
    Constructor.
    """
    Subfield.__init__(self, name)
    return


  def initialize(self, normalizer, spaceDim):
    """
    Initialize subfield metadata.
    """
    from pylith.topology.Field import Field
    self.vectorFieldType = Field.SCALAR
    self.ncomponents = 1
    self.scale = normalizer.pressureScale()
    return


# FACTORIES ////////////////////////////////////////////////////////////

def subfield():
  """
  Factory associated with SubfieldPressure.
  """
  return SubfieldPressure()


# End of file 

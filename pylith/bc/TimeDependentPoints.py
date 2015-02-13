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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/bc/TimeDependentPoints.py
##
## @brief Python abstract base class for managing a boundary condition
## applied to a set of vertices with time dependent paramters.

from TimeDependent import TimeDependent
from bc import TimeDependentPoints as ModuleTimeDependentPoints

def validateDOF(value):
  """
  Validate list of fixed degrees of freedom.
  """
  try:
    size = len(value)
    num = map(int, value)
    for v in num:
      if v < 0:
        raise ValueError
  except:
    raise ValueError, \
          "'bc_dof' must be a zero based list of indices of degrees of " \
          "freedom at a vertex."
  return num
  

# TimeDependentPoints class
class TimeDependentPoints(TimeDependent, ModuleTimeDependentPoints):
  """
  Python abstract base class for managing a boundary condition applied
  to a set of points with time dependent paramters.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  bcDOF = pyre.inventory.list("bc_dof", default=[],
                                 validator=validateDOF)
  bcDOF.meta['tip'] = "Indices of boundary condition DOF " \
      "(0=1st DOF, 1=2nd DOF, etc)."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timedependentpoints"):
    """
    Constructor.
    """
    TimeDependent.__init__(self, name, facility="timedependentpoints")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    TimeDependent._configure(self)

    import numpy
    bcDOF = numpy.array(self.inventory.bcDOF, dtype=numpy.int32)
    ModuleTimeDependentPoints.bcDOF(self, bcDOF)
    return


# End of file 

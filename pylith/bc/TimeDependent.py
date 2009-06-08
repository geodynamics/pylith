#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
#
# ----------------------------------------------------------------------
#

## @file pylith/bc/TimeDependent.py
##
## @brief Python abstract base class for managing a boundary condition
## with time dependent paramters.
##
## This implementation of a time dependent boundary condition creates
## a parameter with the form: 
##   u(x) = u0(x) + v0(x)*(t-tv(x)) + du(x)*a(t-tu(x))
## where
## u0(x) is the initial value given by a spatial database,
## v0(x) is the rate of change given by a spatial database,
## tv(x) is the start time for the rate of change given by a spatial database,
## du(x) is the change in value given by a spatial database,
## tu(x) is the start time for the change in value given by a spatial datasse,
## a(t) is a time history given by a temporal database.

from pyre.components.Component import Component
from bc import TimeDependent as ModuleTimeDependent

from pylith.utils.NullComponent import NullComponent


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
          "'fixed_dof' must be a zero based list of indices of fixed " \
          "degrees of freedom."
  return num
  

# TimeDependent class
class TimeDependent(Component, ModuleTimeDependent):
  """
  Python abstract base class for managing a boundary condition with
  time dependent paramters.

  This implementation of a time dependent boundary condition creates a
  parameter with the form:

  u(x) = u0(x) + v0(x)*(t-tv(x)) + du(x)*a(t-tu(x))

  where
  u0(x) is the initial value given by a spatial database,
  v0(x) is the rate of change given by a spatial database,
  tv(x) is the start time for the rate of change given by a spatial database,
  du(x) is the change in value given by a spatial database,
  tu(x) is the start time for the change in value given by a spatial datasse,
  a(t) is a time history given by a temporal database.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  bcDOF = pyre.inventory.list("fixed_dof", default=[],
                                 validator=validateDOF)
  bcDOF.meta['tip'] = "Indices of boundary condition DOF " \
      "(0=1st DOF, 1=2nd DOF, etc)."

  from spatialdata.spatialdb.SimpleDB import SimpleDB
  dbInitial = pyre.inventory.facility("db_initial", factory=SimpleDB, 
                                      family="spatial_database")
  dbInitial.meta['tip'] = "Database with initial values."
    
  dbRate = pyre.inventory.facility("db_rate", factory=NullComponent, 
                                      family="spatial_database")
  dbRate.meta['tip'] = "Database with rate of change values."
    
  dbChange = pyre.inventory.facility("db_change", factory=NullComponent, 
                                      family="spatial_database")
  dbChange.meta['tip'] = "Database with temporal change in values."
    
  thChange = pyre.inventory.facility("th_change", factory=NullComponent, 
                                      family="temporal_database")
  thChange.meta['tip'] = "Database with time history."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="timedependent"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="time_dependent")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)

    import numpy
    bcDOF = numpy.array(self.inventory.bcDOF, dtype=numpy.int32)
    ModuleTimeDependent.bcDOF(self, bcDOF)

    if isinstance(self.inventory.dbChange, NullComponent):
      if not isinstance(self.inventory.thChange, NullComponent):
        raise ValueError("Cannot provide a time history temporal database "
                         "without a change in value spatial database "
                         "for time dependent boundary condition '%s'." % \
                             self.label)
    
    if not isinstance(self.inventory.dbInitial, NullComponent):
      ModuleTimeDependent.dbInitial(self, self.inventory.dbInitial)
    if not isinstance(self.inventory.dbRate, NullComponent):
      ModuleTimeDependent.dbRate(self, self.inventory.dbRate)
    if not isinstance(self.inventory.dbChange, NullComponent):
      ModuleTimeDependent.dbChange(self, self.inventory.dbChange)
    if not isinstance(self.inventory.thChange, NullComponent):
      ModuleTimeDependent.thChange(self, self.inventory.thChange)
    return


# End of file 

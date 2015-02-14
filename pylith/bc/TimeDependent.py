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

from pylith.utils.PetscComponent import PetscComponent
from bc import TimeDependent as ModuleTimeDependent

from pylith.utils.NullComponent import NullComponent


# TimeDependent class
class TimeDependent(PetscComponent, ModuleTimeDependent):
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

  dbInitial = pyre.inventory.facility("db_initial", factory=NullComponent, 
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
    PetscComponent.__init__(self, name, facility="time_dependent")
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    PetscComponent._configure(self)

    import numpy
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
      ModuleTimeDependent.dbTimeHistory(self, self.inventory.thChange)
    return


# End of file 

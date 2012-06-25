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
# Copyright (c) 2010-2012 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/faults/TractPerturbation.py
##

## @brief Python object for managing parameters for a kinematic
## earthquake sources.
##
## TractPerturbation is responsible for providing the value of
## specified traction at time t over a fault surface.
##
## Factory: eq_kinematic_src

from pylith.bc.TimeDependent import TimeDependent
from faults import TractPerturbation as ModuleTractPerturbation

from pylith.utils.NullComponent import NullComponent

# TractPerturbation class
class TractPerturbation(TimeDependent, ModuleTractPerturbation):
  """
  Python object for managing specified tractions on a fault surface.

  Factory: traction_perturbation
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="tractperturbation"):
    """
    Constructor.
    """
    TimeDependent.__init__(self, name)
    self._createModuleObj()
    self._loggingPrefix = "TrPt "
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    fields = []
    if not isinstance(self.inventory.dbInitial, NullComponent):
      fields += ["traction_initial_value"]
    if not isinstance(self.inventory.dbRate, NullComponent):
      fields += ["traction_rate_of_change", "traction_rate_start_time"]
    if not isinstance(self.inventory.dbChange, NullComponent):
      fields += ["traction_change_in_value", "traction_change_start_time"]

    self.availableFields = \
        {'vertex': {'info': fields,
                    'data': []},
         'cell': {'info': [],
                  'data': []}}
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    ModuleTractPerturbation.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def traction_perturbation():
  """
  Factory associated with TractPerturbation.
  """
  return TractPerturbation()


# End of file 

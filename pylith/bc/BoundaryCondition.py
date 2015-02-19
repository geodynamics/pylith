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

## @file pylith/bc/BoundaryCondition.py
##
## @brief Python abstract base class for managing a boundary condition.
##
## This implementation of a boundary condition applies to a single
## face of an domain and associates both a quadrature scheme with a
## physical boundary condition. Thus, applying different quadrature
## schemes along a face with the same physical boundary condition
## requires two "bc", which can use the same database.
##
## Factory: boundary_condition

from pylith.utils.PetscComponent import PetscComponent
from bc import BoundaryCondition as ModuleBoundaryCondition

# Validator for label
def validateLabel(value):
  """
  Validate label for group/nodeset/pset.
  """
  if 0 == len(value):
    raise ValueError("Label for group/nodeset/pset in mesh not specified.")
  return value


# Validator for direction
def validateDir(value):
  """
  Validate direction.
  """
  msg = "Direction must be a 3 component vector (list)."
  if not isinstance(value, list):
    raise ValueError(msg)
  if 3 != len(value):
    raise ValueError(msg)
  try:
    nums = map(float, value)
  except:
    raise ValueError(msg)
  return value


# BoundaryCondition class
class BoundaryCondition(PetscComponent, ModuleBoundaryCondition):
  """
  Python abstract base class for managing a boundary condition.

  This implementation of a boundary condition applies to a single
  face of an domain.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing BoundaryCondition facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BoundaryCondition facilities and properties.
    ##
    ## \b Properties
    ## @li \b label Label identifier for boundary.
    ##
    ## \b Facilities

    import pyre.inventory

    label = pyre.inventory.str("label", default="", validator=validateLabel)
    label.meta['tip'] = "Label identifier for boundary."

    upDir = pyre.inventory.list("up_dir", default=[0, 0, 1],
		                validator=validateDir)
    upDir.meta['tip'] = "Direction perpendicular to horizontal " \
		        "tangent direction that is not collinear " \
			"with normal direction."

    from pylith.perf.MemoryLogger import MemoryLogger
    perfLogger = pyre.inventory.facility("perf_logger", family="perf_logger",
                                         factory=MemoryLogger)
    perfLogger.meta['tip'] = "Performance and memory logging."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="boundarycondition"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="boundary_condition")
    self._createModuleObj()
    return


  def preinitialize(self, mesh):
    """
    Setup boundary condition.
    """
    import weakref
    self.mesh = weakref.ref(mesh)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize boundary condition.
    """
    ModuleBoundaryCondition.initialize(self, self.mesh(), self.upDir)
    return


  def finalize(self):
    """
    Cleanup.
    """
    self._modelMemoryUse()
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    try:
      PetscComponent._configure(self)
      ModuleBoundaryCondition.label(self, self.inventory.label)
      self.upDir = map(float, self.inventory.upDir)
      self.perfLogger = self.inventory.perfLogger
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring boundary condition "
                       "(%s):\n%s" % (aliases, err.message))
                         
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    raise NotImplementedError, \
          "Please implement _createModuleObj() in derived class."


  def _modelMemoryUse(self):
    """
    Model memory allocation.
    """
    raise NotImplementedError, \
          "Please implement _modelModelUse() in derived class."
    return


# End of file 

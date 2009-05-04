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

from pyre.components.Component import Component
from bc import BoundaryCondition as ModuleBoundaryCondition

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
class BoundaryCondition(Component, ModuleBoundaryCondition):
  """
  Python abstract base class for managing a boundary condition.

  This implementation of a boundary condition applies to a single
  face of an domain.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
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
    ## @li \b db Database of boundary condition parameters

    import pyre.inventory

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Label identifier for boundary."

    upDir = pyre.inventory.list("up_dir", default=[0, 0, 1],
		                validator=validateDir)
    upDir.meta['tip'] = "Direction perpendicular to horizontal " \
		        "tangent direction that is not collinear " \
			"with normal direction."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    db = pyre.inventory.facility("db", factory=SimpleDB, 
                                 family="spatial_database")
    db.meta['tip'] = "Database of boundary condition parameters."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="boundarycondition"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="boundary_condition")
    self._createModuleObj()
    return


  def preinitialize(self, mesh):
    """
    Setup boundary condition.
    """
    self.label(self.inventory.label)
    self.db(self.inventory.db)
    self.mesh = mesh
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize boundary condition.
    """
    ModuleBoundaryCondition.initialize(self, self.mesh, self.upDir)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Component._configure(self)
    self.upDir = map(float, self.inventory.upDir)
    return


  def _createModuleObj(self):
    """
    Call constructor for module object for access to C++ object.
    """
    raise NotImplementedError, \
          "Please implement _createModuleObj() in derived class."


# End of file 

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

## @file pylith/bc/DirichletBC.py
##
## @brief Python object for managing a Dirichlet (prescribed
## displacements) boundary condition with a set of points.
##
## Factory: boundary_condition

from BoundaryCondition import BoundaryCondition
from pylith.feassemble.Constraint import Constraint
from bc import DirichletBC as ModuleDirichletBC

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
  

# DirichletBC class
class DirichletBC(BoundaryCondition, Constraint, ModuleDirichletBC):
  """
  Python object for managing a DirichletBC (prescribed displacements)
  boundary condition.

  Factory: boundary_condition
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(BoundaryCondition.Inventory):
    """
    Python object for managing BoundaryCondition facilities and properties.
    """
    
    ## @class Inventory
    ## Python object for managing BoundaryCondition facilities and properties.
    ##
    ## \b Properties
    ## @li \b fixed_dof Indices of fixed DOF (0=1st DOF, 1=2nd DOF, etc).
    ## @li \b reference_t Reference time for rate of change of values.
    ##
    ## \b Facilities
    ## @li \b initial_db Database of parameters for initial values.
    ## @li \b rate_db Database of parameters for rate of change of values.

    import pyre.inventory

    fixedDOF = pyre.inventory.list("fixed_dof", default=[],
                                   validator=validateDOF)
    fixedDOF.meta['tip'] = "Indices of fixed DOF (0=1st DOF, 1=2nd DOF, etc)."

    from pyre.units.time import s
    tRef = pyre.inventory.dimensional("reference_t", default=0.0*s)
    tRef.meta['tip'] = "Reference time for rate of change of values."

    from FixedDOFDB import FixedDOFDB
    db = pyre.inventory.facility("db", factory=FixedDOFDB,
                                 family="spatial_database")
    db.meta['tip'] = "Database of parameters for initial values."

    dbRate = pyre.inventory.facility("rate_db", factory=FixedDOFDB,
                                 family="spatial_database")
    dbRate.meta['tip'] = "Database of parameters for rate of change of values."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dirichletpoints"):
    """
    Constructor.
    """
    BoundaryCondition.__init__(self, name)
    Constraint.__init__(self)
    self._loggingPrefix = "DiBC "
    return


  def preinitialize(self, mesh):
    """
    Do pre-initialization setup.
    """
    BoundaryCondition.preinitialize(self, mesh)
    Constraint.preinitialize(self, mesh)
    self.dbRate(self.inventory.dbRate)
    self.fixedDOF(self.inventory.fixedDOF)
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    BoundaryCondition.verifyConfiguration(self)
    Constraint.verifyConfiguration(self)

    self._logger.eventEnd(logEvent)
    return


  def initialize(self, totalTime, numTimeSteps, normalizer):
    """
    Initialize DirichletBC boundary condition.
    """
    logEvent = "%sinit" % self._loggingPrefix
    self._logger.eventBegin(logEvent)

    timeScale = normalizer.timeScale()
    tRef = normalizer.nondimensionalize(self.inventory.tRef, timeScale)
    self.referenceTime(tRef)
    
    self.normalizer(normalizer)

    BoundaryCondition.initialize(self, totalTime, numTimeSteps, normalizer)

    self._logger.eventEnd(logEvent)    
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    BoundaryCondition._configure(self)
    return


  def _createModuleObj(self):
    """
    Create handle to corresponding C++ object.
    """
    if None == self.this:
      ModuleDirichletBC.__init__(self)
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def boundary_condition():
  """
  Factory associated with DirichletBC.
  """
  return DirichletBC()

  
# End of file 

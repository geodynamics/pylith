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

## @file pylith/problems/Problem.py
##
## @brief Python abstract base class for crustal dynamics problems.
##
## Factory: problem.

from pyre.components.Component import Component

# Problem class
class Problem(Component):
  """
  Python abstract base class for crustal dynamics problems.

  Factory: problem.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """
    Python object for managing Problem facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Problem facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b materials Materials in problem.
    ## @li \b bc Boundary conditions.
    ## @li \b interfaces Interior surfaces with constraints or
    ##   constitutive models.

    import pyre.inventory
    from pylith.utils.ObjectBin import ObjectBin

    from pylith.materials.Homogeneous import Homogeneous
    materials = pyre.inventory.facility("materials", family="object_bin",
                                        factory=Homogeneous)
    materials.meta['tip'] = "Materials in problem."

    bc = pyre.inventory.facility("bc", family="object_bin", factory=ObjectBin)
    bc.meta['tip'] = "Boundary conditions."

    interfaces = pyre.inventory.facility("interfaces", family="object_bin",
                                         factory=ObjectBin)
    interfaces.meta['tip'] = "Interior surfaces with constraints or " \
                             "constitutive models."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="problem"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="problem")
    self.mesh = None
    return


  def preinitialize(self, mesh):
    """
    Setup integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    raise NotImplementedError, "initialize() not implemented."
    return


  def verifyConfiguration(self):
    """
    Verify compatibility of configuration.
    """
    logEvent = "%sverify" % self._loggingPrefix

    self._logger.eventBegin(logEvent)
    for material in self.materials.bin:
      material.verifyConfiguration()
    for bc in self.bc.bin:
      bc.verifyConfiguration()
    for interface in self.interfaces.bin:
      interface.verifyConfiguration()
    self._logger.eventEnd(logEvent)
    return
  

  def initialize(self):
    """
    Initialize integrators for each element family (material/quadrature,
    bc/quadrature, etc.).
    """
    raise NotImplementedError, "initialize() not implemented."
    return


  def run(self, app):
    """
    Solve the problem.
    """
    raise NotImplementedError, "run() not implemented."
    return


  def finalize(self, mesh):
    """
    Cleanup.
    """
    raise NotImplementedError, "finalize() not implemented."
    return


  def checkpoint(self):
    """
    Save problem state for restart.
    """
    raise NotImplementedError, "checkpoint() not implemented."
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.materials = self.inventory.materials
    self.bc = self.inventory.bc
    self.interfaces = self.inventory.interfaces
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    if not "_loggingPrefix" in dir(self):
      self._loggingPrefix = ""

    from pylith.utils.EventLogger import EventLogger
    logger = EventLogger()
    logger.setClassName("Problem")
    logger.initialize()

    events = ["preinit",
              "verify",
              "init",
              "run",
              "finalize"]
    for event in events:
      logger.registerEvent("%s%s" % (self._loggingPrefix, event))

    self._logger = logger
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def problem():
  """
  Factory associated with Problem.
  """
  return Problem()


# End of file 

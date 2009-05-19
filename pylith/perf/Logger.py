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

## @file pylith/perf/Logger.py
##
## @brief Python abstract base class for performance and memory logging.
##
## Factory: perf_logger.

from pylith.utils.PetscComponent import PetscComponent

# Logger class
class Logger(PetscComponent):
  """
  Python abstract base class for performance and memory logging.

  Factory: perf_logger.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing Problem facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Problem facilities and properties.
    ##
    ## \b Properties
    ## @li \b verbose Should information be printed to the screen.

    import pyre.inventory

    verbose = pyre.inventory.bool("verbose", default=False)
    verbose.meta['tip'] = "Print information to the screen."

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="perf_logger"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="perf_logger")
    return

  def join(self, logger):
    """
    Incorporate information from another logger.
    """
    raise NotImplementedError, "join() not implemented."
    return

  def logMesh(self, stage, mesh):
    """
    Read mesh parameters to determine memory from our model.
    """
    raise NotImplementedError, "logMesh() not implemented."
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.verbose = self.inventory.verbose
    return


  def _setupLogging(self):
    """
    Setup event logging.
    """
    return
  

# FACTORIES ////////////////////////////////////////////////////////////

def perf_logger():
  """
  Factory associated with Logger.
  """
  return Logger()


# End of file 

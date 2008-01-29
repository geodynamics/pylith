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

## @file pylith/bc/FixedDOFDB.py
##
## @brief Python object for spatial database with uniform zero values
## for degrees of freedom.
##
## Factory: spatial_database

from spatialdata.spatialdb.UniformDB import UniformDB

# FixedDOFDB class
class FixedDOFDB(UniformDB):
  """
  Python object for spatial database with uniform zero values for
  degrees of freedom.

  Factory: spatial_database
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(UniformDB.Inventory):
    """
    Python object for managing FixedDOFDB facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing FixedDOFDB facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li none

    import pyre.inventory

    values = ["dof-0", "dof-1", "dof-2"]
    data = [0.0, 0.0, 0.0]


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fixeddofdb"):
    """
    Constructor.
    """
    UniformDB.__init__(self, name)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based on inventory.
    """
    UniformDB._configure(self)
    return


# FACTORIES ////////////////////////////////////////////////////////////

def spatial_database():
  """
  Factory associated with FixedDOFDB.
  """
  return FixedDOFDB()


# End of file 

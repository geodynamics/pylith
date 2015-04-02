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

## @file pylith/bc/ZeroDispDB.py
##
## @brief Python object for spatial database with uniform zero
## displacements for degrees of freedom.
##
## Factory: spatial_database

from spatialdata.spatialdb.UniformDB import UniformDB

# ZeroDispDB class
class ZeroDispDB(UniformDB):
  """
  Python object for spatial database with uniform zero displacements
  for degrees of freedom.

  Factory: spatial_database
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(UniformDB.Inventory):
    """
    Python object for managing ZeroDispDB facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ZeroDispDB facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li none

    import pyre.inventory

    from pyre.units.length import m
    values = ["displacement-x", "displacement-y", "displacement-z"]
    data = [0.0*m, 0.0*m, 0.0*m]

    label = pyre.inventory.str("label", default="Zero displacement BC.")


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="zerodispdb"):
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
  Factory associated with ZeroDispDB.
  """
  return ZeroDispDB()


# End of file 

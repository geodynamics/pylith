#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
# Romain Jolivet, Caltech
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

## @file pylith/faults/DKSelector.py
##
## @brief Python object for the dynamic-kinematic selector
##

from faults import DKSelector as ModuleDKSelector

# DKSelector class
class DKSelector(ModuleDKSelector):
  """
  Python object for a dynamic-kinematic selector

  Inventory

  \b Properties
  @li None
  
  \b Facilities
  @li \b dynamic_kinematic_selector Spatial database of dk selector

  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  from spatialdata.spatialdb.SimpleDB import SimpleDB
  
  dbDKSel = pyre.inventory.facility("dynamic_kinematic_selector", family="spatial_database",
                                       factory=SimpleDB)
  dbDKSel.meta['tip'] = "Spatial database of dynamic-kinematic selector"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="dkselector"):
    """
    Constructor.
    """
    ModuleDKSelector.__init__(self)
    self._loggingPrefix = "DKsel "
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    ModuleDKSelector.dbDKSel(self, self.inventory.dbDKSel)
    return

# End of file 

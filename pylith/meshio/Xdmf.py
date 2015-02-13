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

## @file pyre/meshio/Xdmf.py
##
## @brief Python class for Xdmf metadata file associated with an HDF5 file.
##
## Factory: xdmf

from pylith.utils.PetscComponent import PetscComponent
from meshio import Xdmf as ModuleXdmf

# Xdmf class
class Xdmf(PetscComponent, ModuleXdmf):
  """
  Python class for Xdmf metadata file associated with an HDF5 file.

  \b Properties
  @li \b filename Name of HDF5 file.
  
  \b Facilities
  @li None
  """

  # INVENTORY //////////////////////////////////////////////////////////

  import pyre.inventory

  filename = pyre.inventory.str("filename", default="output.h5")
  filename.meta['tip'] = "Name of HDF5 file."
  

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="xdmf"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="xdmf")
    self._loggingPrefix = "Xdmf "

    ModuleXdmf.__init__(self)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.filename = self.inventory.filename
    return


# End of file 

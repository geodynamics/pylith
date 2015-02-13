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


## @file pylith/topology/JacobianViewer.py
##
## @brief Python object for writing system Jacobian to file.
##
## Factory: jacobian_viewer

from pylith.utils.PetscComponent import PetscComponent

# JacobianViewer class
class JacobianViewer(PetscComponent):
  """
  Python abstract base class for formulations of solving equations.

  In general, we use some explicit or implicit formulation of the PDEs
  to create a linear form, [A]{u}={b} that we can solve.

  Factory: pde_formulation.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """
    Python object for managing JacobianViewer facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing JacobianViewer facilities and properties.
    ##
    ## \b Properties
    ## @li \b filename Filename for Jacobian matrix.
    ## @li \b time_format C style format string for time stamp in filename.
    ## @li \b time_constant Value used to normalize time stamp in filename.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    filename = pyre.inventory.str("filename", default="jacobian.mat")
    filename.meta['tip'] = "Filename for Jacobian matrix."

    timeFormat = pyre.inventory.str("time_format", default="%f")
    timeFormat.meta['tip'] = "C style format string for time stamp in filename."

    from pyre.units.time import second
    timeConstant = pyre.inventory.dimensional("time_constant",
                                              default=1.0*second,
                              validator=pyre.inventory.greater(0.0*second))
    timeConstant.meta['tip'] = \
        "Values used to normalize time stamp in filename."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="formulation"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="jacobian_viewer")
    return


  def view(self, jacobian, t, comm):
    """
    Write Jacobian to binary file.
    """
    jacobian.write(self._filenameStamp(t), comm)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.filename = self.inventory.filename
    self.timeFormat = self.inventory.timeFormat
    self.timeConstant = self.inventory.timeConstant
    return


  def _filenameStamp(self, t):
    """
    Create filename by extracting basename and adding a time stamp.
    """
    timeStamp = self.timeFormat % (t/self.timeConstant.value)
    basename = self.filename
    if basename.endswith(".mat"):
      basename = basename[0:len(basename)-4]
    filename = basename + "_t" + timeStamp + ".mat"
    return filename

      
# FACTORIES ////////////////////////////////////////////////////////////

def jacobian_viewer():
  """
  Factory associated with JacobianViewer.
  """
  return JacobianViewer()


# End of file 

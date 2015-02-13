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

## @file pylith/feassemble/quadrature/Qudrature.py
##
## @brief Python abstract base class for integrating over
## finite-elements using quadrature.

from pylith.utils.PetscComponent import PetscComponent

# ----------------------------------------------------------------------
# QuadratureBase class
class QuadratureBase(PetscComponent):
  """
  Python abstract base class for integrating over finite-elements
  using quadrature.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(PetscComponent.Inventory):
    """Python object for managing Quadrature facilities and properties."""

    ## @class Inventory
    ## Python object for managing Quadrature facilities and properties.
    ##
    ## \b Properties
    ## @li \b min_jacobian Minimum allowable determinant of Jacobian.
    ## @li \b check_conditoning Check element matrices for 
    ##   ill-conditioning.
    ##
    ## \b Facilities
    ## @li \b cell Reference cell with basis functions and quadrature rules

    import pyre.inventory

    minJacobian = pyre.inventory.float("min_jacobian", default=1.0e-06)
    minJacobian.meta['tip'] = "Minimum allowable determinant of Jacobian."

    checkConditioning = pyre.inventory.bool("check_conditioning",
                                            default=False)
    checkConditioning.meta['tip'] = \
        "Check element matrices for ill-conditioning."

    from pylith.feassemble.FIATSimplex import FIATSimplex
    cell = pyre.inventory.facility("cell", family="reference_cell",
                                   factory=FIATSimplex)
    cell.meta['tip'] = "Reference cell with basis fns and quadrature rules."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature"):
    """
    Constructor.
    """
    PetscComponent.__init__(self, name, facility="quadrature")
    return


  def preinitialize(self, spaceDim):
    """
    Setup quadrature object.
    """
    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    if 0 == comm.rank:
      self._info.log("Initializing reference cell.")
    cell = self.cell
    cell.initialize(spaceDim)

    if 0 == comm.rank:
      self._info.log("Initializing C++ quadrature.")
    self._initialize(cell)
    self.refGeometry(cell.geometry)
    return


  def initialize(self):
    """
    Initialize quadrature object.
    """
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    PetscComponent._configure(self)
    self.minJacobian(self.inventory.minJacobian)
    self.checkConditioning(self.inventory.checkConditioning)
    self.cell = self.inventory.cell
    return


# ----------------------------------------------------------------------
from feassemble import Quadrature as ModuleQuadrature

# Quadrature class
class Quadrature(QuadratureBase, ModuleQuadrature):
  """
  Python object for integrating over finite-elements using quadrature.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature"):
    """
    Constructor.
    """
    QuadratureBase.__init__(self, name)
    ModuleQuadrature.__init__(self)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _initialize(self, cell):
    """
    Initialize C++ quadrature object.
    """
    import numpy
    from pylith.utils.utils import sizeofPylithScalar
    size = sizeofPylithScalar()
    if 8 == size:
      ModuleQuadrature.initialize(self, cell.basis,
                                  cell.basisDeriv,
                                  cell.quadPts,
                                  cell.quadWts,
                              cell.geometry.spaceDim())
    elif 4 == size:
      ModuleQuadrature.initialize(self, numpy.float32(cell.basis),
                                  numpy.float32(cell.basisDeriv),
                                  numpy.float32(cell.quadPts),
                                  numpy.float32(cell.quadWts),
        cell.geometry.spaceDim())
    else:
        raise ValueError("Unknown size for PylithScalar")
    return


# End of file 

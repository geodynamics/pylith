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

## @file pylith/feassemble/Qudrature.py

## @brief Python abstract base class for integrating over
## finite-elements using quadrature.

from pyre.components.Component import Component

# ----------------------------------------------------------------------
# Quadrature class
class Quadrature(Component):
  """
  Python abstract base class for integrating over finite-elements
  using quadrature.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Quadrature facilities and properties."""

    ## @class Inventory
    ## Python object for managing Quadrature facilities and properties.
    ##
    ## \b Properties
    ## @li \b min_jacobian Minimum allowable determinant of Jacobian.
    ##
    ## \b Facilities
    ## @li \b cell Reference cell with basis functions and quadrature rules

    import pyre.inventory

    minJacobian = pyre.inventory.float("min_jacobian", default=1.0e-06)
    minJacobian.meta['tip'] = "Minimum allowable determinant of Jacobian."

    from FIATSimplex import FIATSimplex
    cell = pyre.inventory.facility("cell", factory=FIATSimplex)
    cell.meta['tip'] = "Reference cell with basis fns and quadrature rules."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="quadrature"):
    """
    Constructor.
    """
    Component.__init__(self, name, facility="quadrature")
    self.cppHandle = None
    self.spaceDim = None
    return

  def initialize(self):
    """
    Initialize C++ quadrature object.
    """
    self.cppHandle.minJacobian = self.minJacobian

    c = self.cell
    c.initialize()
    self.cppHandle.initialize(c.basis, c.basisDeriv,
                              c.quadrature.get_points(),
                              c.quadrature.get_weights(),
                              c.cellDim, c.numCorners, c.numQuadPts,
                              self.spaceDim)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    Component._configure(self)
    self.minJacobian = self.inventory.minJacobian
    self.cell = self.inventory.cell
    return


# ----------------------------------------------------------------------
# Quadrature1D class
class Quadrature1D(Quadrature):
  """
  Python object for integrating over 1-D finite-elements in a 1-D
  domain using quadrature.
  """

  def __init__(self, name="quadrature1d"):
    """
    Constructor.
    """
    Quadrature.__init(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature1D()
    self.spaceDim = 1
    return


# ----------------------------------------------------------------------
# Quadrature1Din2D class
class Quadrature1Din2D(Quadrature):
  """
  Python object for integrating over 1-D finite-elements in a 2-D
  domain using quadrature.
  """

  def __init__(self, name="quadrature1din2d"):
    """
    Constructor.
    """
    Quadrature.__init(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature1Din2D()
    self.spaceDim = 2
    return


# ----------------------------------------------------------------------
# Quadrature1Din3D class
class Quadrature1Din3D(Quadrature):
  """
  Python object for integrating over 1-D finite-elements in a 3-D
  domain using quadrature.
  """

  def __init__(self, name="quadrature1din3d"):
    """
    Constructor.
    """
    Quadrature.__init(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature1Din3D()
    self.spaceDim = 3
    return


# ----------------------------------------------------------------------
# Quadrature2D class
class Quadrature2D(Quadrature):
  """
  Python object for integrating over 2-D finite-elements in a 2-D
  domain using quadrature.
  """

  def __init__(self, name="quadrature2d"):
    """
    Constructor.
    """
    Quadrature.__init(self, name)
    import pylith.feassemble.feassemble as bindings
    self.cppHandle = bindings.Quadrature@D()
    self.spaceDim = 2
    return


# End of file 

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

## @file pylith/feassemble/FIATSimplex.py
##
## @brief Python object for managing basis functions and quadrature
## rules of a simplex reference finite-element cell using FIAT.
##
## Factory: reference_cell.

from FIATCell import FIATCell

import numpy

def validateShape(shape):
  name = shape.lower()
  if not ("tetrahedron" == name or 
          "triangle" == name or 
          "line" == name):
    raise ValueError("Unknown shape '%s' for reference finite-element " \
                     "cell.\n" \
                     "Known shapes: 'tetrahedron', 'triangle', 'line'" % \
                     name)
  return name

# FIATSimplex class
class FIATSimplex(FIATCell):
  """
  Python object for managing basis functions and quadrature rules of a
  simplex reference finite-element cell using FIAT.

  Factory: reference_cell.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(FIATCell.Inventory):
    """Python object for managing FIATSimplex facilities and properties."""

    ## @class Inventory
    ## Python object for managing FIATSimplex facilities and properties.
    ##
    ## \b Properties
    ## @li \b shape Shape of finite-element cell
    ## @li \b degree Degree of finite-element cell 
    ## @li \b quad_order Order of quadrature rule
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    shape = pyre.inventory.str("shape", default="tetrahedron",
                               validator=validateShape)
    shape.meta['tip'] = "Shape of finite-element cell."

    degree = pyre.inventory.int("degree", default=1)
    degree.meta['tip'] = "Degree of finite-element cell."

    order = pyre.inventory.int("quad_order", default=3)
    order.meta['tip'] = "Order of quadrature rule."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fiatsimplex"):
    """
    Constructor.
    """
    FIATCell.__init__(self, name)
    return

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    FIATCell._configure(self)
    self.shape = self.inventory.shape
    self.degree = self.inventory.degree
    self.order = self.inventory.order
    return


  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    import FIAT.quadrature
    return FIAT.quadrature.make_quadrature(self._getShape(), self.order)


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    from FIAT.Lagrange import Lagrange
    return Lagrange(self._getShape(), self.degree).function_space()


  def _setupVertices(self):
    """
    Setup vertices for reference cell.
    """
    import FIAT.shapes
    return FIAT.shapes.vertices[self._getShape()].values()


  def _getShape(self):
    """
    Parse string into FIAT shape.
    """
    import FIAT.shapes
    name = self.shape.lower()
    if "tetrahedron" == name:
      shape = FIAT.shapes.TETRAHEDRON
    elif "triangle" == name:
      shape = FIAT.shapes.TRIANGLE
    elif "line" == name:
      shape = FIAT.shapes.LINE
    else:
      raise ValueError("Unknown shape '%s' for reference finite-element " \
                       "cell.\n" \
                       "Known shapes: 'tetrahedron', 'triangle', 'line'" % \
                       name)
    return shape


# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATSimplex.
  """
  return FIATSimplex()


# End of file 

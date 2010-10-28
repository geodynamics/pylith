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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/feassemble/FIATSimplex.py
##
## @brief Python object for managing basis functions and quadrature
## rules of a simplex reference finite-element cell using FIAT.
##
## Factory: reference_cell.

from ReferenceCell import ReferenceCell

import numpy

def validateShape(shape):
  name = shape.lower()
  if not ("tetrahedron" == name or 
          "triangle" == name or 
          "line" == name or
          "point" == name):
    raise ValueError("Unknown shape '%s' for reference finite-element " \
                     "cell.\n" \
                     "Known shapes: 'tetrahedron', 'triangle', 'line', 'point'" % \
                     name)
  return name

# FIATSimplex class
class FIATSimplex(ReferenceCell):
  """
  Python object for managing basis functions and quadrature rules of a
  simplex reference finite-element cell using FIAT.

  Factory: reference_cell.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ReferenceCell.Inventory):
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

    order = pyre.inventory.int("quad_order", default=-1)
    order.meta['tip'] = "Order of quadrature rule [-1, order = degree]."
    

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fiatsimplex"):
    """
    Constructor.
    """
    ReferenceCell.__init__(self, name)
    return


  def initialize(self, spaceDim):
    """
    Initialize reference finite-element cell.
    """
    self._setupGeometry(spaceDim)

    if "point" == self.shape.lower():
      # Need 0-D quadrature for boundary conditions of 1-D meshes
      self.cellDim = 0
      self.numCorners = 1
      self.numQuadPts = 1
      self.basis = numpy.array([[1.0]])
      self.basisDeriv = numpy.array([[[1.0]]])
      self.quadPts = numpy.array([[0.0]])
      self.quadWts = numpy.array([1.0])
      self.vertices = numpy.array([[0.0]])
    else:
      quadrature = self._setupQuadrature()
      basisFns = self._setupBasisFns()

      # Get coordinates of vertices (dual basis)
      vertices = numpy.array(self._setupVertices(), dtype=numpy.float64)

      # Evaluate basis functions at quadrature points
      quadpts = quadrature.get_points()
      basis = numpy.array(basisFns.tabulate(quadpts)).transpose()

      # Evaluate derivatives of basis functions at quadrature points
      import FIAT.shapes
      dim = FIAT.shapes.dimension(basisFns.base.shape)
      basisDeriv = numpy.array([basisFns.deriv_all(d).tabulate(quadpts) \
                                for d in range(dim)]).transpose()

      self.cellDim = dim
      self.numCorners = len(basisFns)
      self.numQuadPts = len(quadrature.get_weights())

      # Permute from FIAT order to Sieve order
      p = self._permutationFIATToSieve()
      self.vertices = vertices[p,:]
      self.basis = numpy.reshape(basis[:,p].flatten(), basis.shape)
      self.basisDeriv = numpy.reshape(basisDeriv[:,p,:].flatten(), 
                                      basisDeriv.shape)

      # No permutation in order of quadrature points
      self.quadPts = numpy.array(quadrature.get_points())
      self.quadWts = numpy.array(quadrature.get_weights())


    self._info.line("Cell geometry: ")
    self._info.line(self.geometry)
    self._info.line("Vertices: ")
    self._info.line(self.vertices)
    self._info.line("Quad pts:")
    self._info.line(self.quadPts)
    self._info.line("Quad wts:")
    self._info.line(self.quadWts)
    self._info.line("Basis fns @ quad pts ):")
    self._info.line(self.basis)
    self._info.line("Basis fn derivatives @ quad pts:")
    self._info.line(self.basisDeriv)
    self._info.log()    
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members based using inventory.
    """
    ReferenceCell._configure(self)
    self.shape = self.inventory.shape
    self.degree = self.inventory.degree
    self.order = self.inventory.order
    if self.order == -1:
      self.order = self.degree
    return

  
  def _permutationFIATToSieve(self):
    """
    Permute from FIAT basis order to Sieve basis order.

    FIAT: corners, edges, faces

    Sieve: breadth search first (faces, edges, corners)
    """
    import FIAT.shapes

    basis = self.cell.function_space()
    dim = FIAT.shapes.dimension(self._getShape())
    ids = self.cell.Udual.entity_ids
    permutation = []
    if dim == 1:
      for vertex in ids[0]:
        permutation.extend(ids[0][vertex])
      for edge in ids[1]:
        permutation.extend(ids[1][edge])
    elif dim == 2:
      for vertex in ids[0]:
        permutation.extend(ids[0][vertex])
      for edge in ids[1]:
        permutation.extend(ids[1][(edge+2)%3])
      for face in ids[2]:
        permutation.extend(ids[2][face])
    elif dim == 3:
      for vertex in ids[0]:
        permutation.extend(ids[0][vertex])
      for edge in [2, 0, 1, 3]:
        permutation.extend(ids[1][edge])
      for edge in [4, 5]:
        if len(ids[1][edge]) > 0:
          permutation.extend(ids[1][edge][::-1])
      for face in [3, 2, 0, 1]:
        permutation.extend(ids[2][face])
      for volume in ids[3]:
        permutation.extend(ids[3][volume])
    else:
      raise ValueError("Unknown dimension '%d'." % dim)
    return permutation


  def _setupGeometry(self, spaceDim):
    """
    Setup reference cell geometry object.
    """
    import CellGeometry

    self.geometry = None
    name = self.shape.lower()
    if "tetrahedron" == name:
      if 3 == spaceDim:
        self.geometry = CellGeometry.GeometryTet3D()
    elif "triangle" == name:
      if 2 == spaceDim:
        self.geometry = CellGeometry.GeometryTri2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryTri3D()
    elif "line" == name:
      if 1 == spaceDim:
        self.geometry = CellGeometry.GeometryLine1D()
      elif 2 == spaceDim:
        self.geometry = CellGeometry.GeometryLine2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryLine3D()
    elif "point" == name:
      if 1 == spaceDim:
        self.geometry = CellGeometry.GeometryPoint1D()
      elif 2 == spaceDim:
        self.geometry = CellGeometry.GeometryPoint2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryPoint3D()
    if None == self.geometry:
      raise ValueError("Could not set shape '%s' of cell for '%s' in spatial " \
                       "dimension '%s'." % (name, self.name, spaceDim))

    return
  

  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    
    import FIAT.quadrature
    return FIAT.quadrature.make_quadrature_by_degree(self._getShape(),
                                                     self.order)


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    from FIAT.Lagrange import Lagrange
    self.cell = Lagrange(self._getShape(), self.degree)
    return self.cell.function_space() 


  def _setupVertices(self):
    """
    Setup evaluation functional points for reference cell.
    """
    return self.cell.Udual.pts


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
    elif "point" == name:
      shape = None
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

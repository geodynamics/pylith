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

## @file pylith/feassemble/FIATSimplex.py
##
## @brief Python object for managing basis functions and quadrature
## rules of a simplex reference finite-element cell using FIAT.
##
## Factory: reference_cell.

from ReferenceCell import ReferenceCell

import numpy

def validateDimension(dim):
  if dim < 0 or dim > 3:
    raise ValueError("Dimension of simplex cell must be 0, 1, 2, or 3.")
  return dim


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
    ## @li \b dimension Dimension of finite-element cell.
    ## @li \b degree Degree of finite-element cell 
    ## @li \b quad_order Order of quadrature rule
    ## @li \b collocate_quad Collocate quadrature points with vertices.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    dimension = pyre.inventory.int("dimension", default=3,
                                   validator=validateDimension)
    dimension.meta['tip'] = "Dimension of finite-element cell."

    degree = pyre.inventory.int("degree", default=1)
    degree.meta['tip'] = "Degree of finite-element cell."

    order = pyre.inventory.int("quad_order", default=-1)
    order.meta['tip'] = "Order of quadrature rule [-1, order = degree]."
    
    collocateQuad = pyre.inventory.bool("collocate_quad", default=False)
    collocateQuad.meta['tip'] = "Collocate quadrature points with vertices."
    

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

    quadrature = self._setupQuadrature()
    basisFns = self._setupBasisFns()

    # Get coordinates of vertices (dual basis)
    vertices = numpy.array(self._setupVertices(), dtype=numpy.float64)

    # Evaluate basis functions at quadrature points
    from FIAT.polynomial_set import mis
    quadpts = quadrature.get_points()
    dim     = basisFns.ref_el.get_spatial_dimension()
    evals   = basisFns.tabulate(quadpts, 1)
    basis   = numpy.array(evals[mis(dim, 0)[0]], dtype=numpy.float64).transpose()

    # Evaluate derivatives of basis functions at quadrature points
    basisDeriv = numpy.array([evals[alpha] for alpha in mis(dim, 1)], dtype=numpy.float64).transpose()

    self.cellDim = dim
    self.numCorners = basisFns.get_num_members()
    self.numQuadPts = len(quadrature.get_weights())

    # Permute from FIAT order to Sieve order
    p = self._permutationFIATToSieve()
    self.vertices = vertices[p,:]
    self.basis = numpy.reshape(basis[:,p].flatten(), basis.shape)
    self.basisDeriv = numpy.reshape(basisDeriv[:,p,:].flatten(), 
                                    basisDeriv.shape)

    # No permutation in order of quadrature points
    self.quadPts = numpy.array(quadrature.get_points(), dtype=numpy.float64)
    self.quadWts = numpy.array(quadrature.get_weights(), dtype=numpy.float64)


    from pylith.mpi.Communicator import mpi_comm_world
    comm = mpi_comm_world()
    if 0 == comm.rank:
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
    try:
      ReferenceCell._configure(self)
      self.cellDim = self.inventory.dimension
      self.degree = self.inventory.degree
      self.order = self.inventory.order
      self.collocateQuad = self.inventory.collocateQuad

      if self.order == -1:
        self.order = self.degree
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring quadrature "
                       "(%s):\n%s" % (aliases, err.message))
    return

  
  def _permutationFIATToSieve(self):
    """
    Permute from FIAT basis order to Sieve basis order.

    FIAT: corners, edges, faces

    Sieve: breadth search first (faces, edges, corners)
    """
    basis = self.cell.get_nodal_basis()
    dim = self._getShape().get_spatial_dimension()
    ids = self.cell.dual.get_entity_ids()
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
      # Flip vertices 0 and 1
      vids = [v for v in ids[0]]
      tmp     = vids[0]
      vids[0] = vids[1]
      vids[1] = tmp
      for vertex in vids:
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
    if self.cellDim == 3:
      if 3 == spaceDim:
        self.geometry = CellGeometry.GeometryTet3D()
    elif self.cellDim == 2:
      if 2 == spaceDim:
        self.geometry = CellGeometry.GeometryTri2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryTri3D()
    elif self.cellDim == 1:
      if 2 == spaceDim:
        self.geometry = CellGeometry.GeometryLine2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryLine3D()
    if None == self.geometry:
      raise ValueError("Could not set geometry of cell with dimension '%d' for '%s' in spatial " \
                       "dimension '%s'." % (self.cellDim, self.name, spaceDim))

    return
  

  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    from FIAT.reference_element import default_simplex
    from FIAT.quadrature import make_quadrature
    from FIATQuadrature import CollocatedQuadratureRule
      
    if not self.collocateQuad:
      q = make_quadrature(self._getShape(), self.order)
    else:
      q = CollocatedQuadratureRule(self._getShape(), self.order)

    return q


  def _setupBasisFns(self):
    """
    Setup basis functions for reference cell.
    """
    from FIAT.lagrange import Lagrange
    self.cell = Lagrange(self._getShape(), self.degree)
    return self.cell.get_nodal_basis() 


  def _setupVertices(self):
    """
    Setup evaluation functional points for reference cell.
    """
    return numpy.array([n.get_point_dict().keys()[0] for n in self.cell.dual.get_nodes()], dtype=numpy.float64)


  def _getShape(self):
    """
    Parse string into FIAT shape.
    """
    from FIAT.reference_element import default_simplex
    if self.cellDim == 3:
      shape = default_simplex(3)
    elif self.cellDim == 2:
      shape = default_simplex(2)
    elif self.cellDim == 1:
      shape = default_simplex(1)
    elif self.cellDim == 0:
      shape = None
    else:
      raise ValueError("Unknown dimension '%d' for reference finite-element " \
                       "cell.\n" % self.cellDim)
    return shape


# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATSimplex.
  """
  return FIATSimplex()


# End of file 

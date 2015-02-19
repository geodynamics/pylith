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

## @file pylith/feassemble/FIATLagrange.py
##
## @brief Python object for managing basis functions and quadrature
## rules of a Lagrange reference finite-element cell using FIAT.
##
## The basis functions are constructed from the tensor product of 1-D
## Lagrance reference cells.
##
## Factory: reference_cell.

from ReferenceCell import ReferenceCell

import numpy

def validateDimension(dim):
  if dim < 0 or dim > 3:
    raise ValueError("Dimension of Lagrange cell must be 0, 1, 2, or 3.")
  return dim


# FIATLagrange class
class FIATLagrange(ReferenceCell):
  """
  Python object for managing basis functions and quadrature rules of a
  Lagrange reference finite-element cell using FIAT.

  Factory: reference_cell.
  """

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(ReferenceCell.Inventory):
    """Python object for managing FIATLagrange facilities and properties."""

    ## @class Inventory
    ## Python object for managing FIATLagrange facilities and properties.
    ##
    ## \b Properties
    ## @li \b dimension Dimension of finite-element cell.
    ## @li \b degree Degree of finite-element cell.
    ## @li \b quad_order Order of quadrature rule.
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
    order.meta['tip'] = "Order of quadrature rule."
    
    collocateQuad = pyre.inventory.bool("collocate_quad", default=False)
    collocateQuad.meta['tip'] = "Collocate quadrature points with vertices."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="fiatlagrange"):
    """
    Constructor.
    """
    ReferenceCell.__init__(self, name)
    return


  def initialize(self, spaceDim):
    """
    Initialize reference finite-element cell from a tensor product of
    1-D Lagrange elements.
    """
    self._setupGeometry(spaceDim)

    if self.cellDim > 0:
      quadrature = self._setupQuadrature()
      element = self._setupElement()
      dim = self.cellDim
    
      # Get coordinates of vertices (dual basis)
      vertices = numpy.array(self._setupVertices(element), dtype=numpy.float64)

      # Evaluate basis functions at quadrature points
      from FIAT.polynomial_set import mis
      quadpts    = numpy.array(quadrature.get_points(), dtype=numpy.float64)
      quadwts    = numpy.array(quadrature.get_weights(), dtype=numpy.float64)
      numQuadPts = len(quadpts)
      evals      = element.get_nodal_basis().tabulate(quadrature.get_points(), 1)
      basis      = numpy.array(evals[mis(1, 0)[0]], dtype=numpy.float64).transpose()
      numBasis   = element.get_nodal_basis().get_num_members()

      # Evaluate derivatives of basis functions at quadrature points
      basisDeriv = numpy.array([evals[alpha] for alpha in mis(1, 1)], dtype=numpy.float64).transpose()

      self.numQuadPts = numQuadPts**dim
      self.numCorners = numBasis**dim

      if dim == 1:
        # Set order of vertices and basis functions.
        # Corners
        vertexOrder = [0, 1]
        # Edges
        p = numpy.arange(2, numBasis, dtype=numpy.int32)
        vertexOrder += zip(p)
        
        self.vertices = numpy.zeros((self.numCorners, dim))
        n = 0
        for p in vertexOrder:
          self.vertices[n][0] = vertices[p]
          n += 1
        if not n == self.numCorners:
          raise RuntimeError('Invalid 1-D vertex ordering: '+str(n)+ \
                             ' should be '+str(self.numCorners))
          
        self.quadPts = numpy.zeros((numQuadPts, dim))
        self.quadWts = numpy.zeros((numQuadPts,))
        self.basis = numpy.zeros((numQuadPts, numBasis))
        self.basisDeriv = numpy.zeros((numQuadPts, numBasis, dim))
        
        # Order of quadrature points doesn't matter
        # Order of basis functions should match vertices for isoparametric
        n = 0
        for p in range(0, numQuadPts):
          self.quadPts[n][0] = quadpts[p]
          self.quadWts[n]    = quadwts[p]
            
          m = 0
          for (bp) in vertexOrder:
            self.basis[n][m] = basis[p][bp]
            self.basisDeriv[n][m][0] = basisDeriv[p][bp][0]
            m += 1
          if not m == self.numCorners:
            raise RuntimeError('Invalid 2-D basis tabulation: '+str(m)+ \
                               ' should be '+str(self.numCorners))
          n += 1
        if not n == self.numQuadPts:
          raise RuntimeError('Invalid 2-D quadrature: '+str(n)+ \
                             ' should be '+str(self.numQuadPts))

      elif dim == 2:
          # Set order of vertices and basis functions.
          # Corners
          vertexOrder = [(0,0), (1,0), (1,1), (0,1)]
          # Edges
          #   Bottom
          p = numpy.arange(2, numBasis, dtype=numpy.int32)
          q = numpy.zeros(numBasis-2, dtype=numpy.int32)
          vertexOrder += zip(p,q)
          #   Right
          p = numpy.ones(numBasis-2, dtype=numpy.int32)
          q = numpy.arange(2, numBasis, dtype=numpy.int32)
          vertexOrder += zip(p,q)
          #   Top
          p = numpy.arange(numBasis-1, 1, step=-1, dtype=numpy.int32)
          q = numpy.ones(numBasis-2, dtype=numpy.int32)
          vertexOrder += zip(p,q)
          #   Left
          p = numpy.zeros(numBasis-2, dtype=numpy.int32)
          q = numpy.arange(numBasis-1, 1, step=-1, dtype=numpy.int32)
          vertexOrder += zip(p,q)
          # Face
          p = numpy.arange(2, numBasis, dtype=numpy.int32)
          q = numpy.arange(2, numBasis, dtype=numpy.int32)
          vertexOrder += zip(p,q)
          
          self.vertices = numpy.zeros((self.numCorners, dim))
          n = 0
          for (p,q) in vertexOrder:
            self.vertices[n][0] = vertices[p]
            self.vertices[n][1] = vertices[q]
            n += 1
          if not n == self.numCorners:
            raise RuntimeError('Invalid 2-D vertex ordering: '+str(n)+ \
                               ' should be '+str(self.numCorners))
        
          self.quadPts = numpy.zeros((numQuadPts*numQuadPts, dim))
          self.quadWts = numpy.zeros((numQuadPts*numQuadPts,))
          self.basis = numpy.zeros((numQuadPts*numQuadPts,
                                         numBasis*numBasis))
          self.basisDeriv = numpy.zeros((numQuadPts*numQuadPts,
                                         numBasis*numBasis, dim))

          # Order of quadrature points doesn't matter
          # Order of basis functions should match vertices for isoparametric
          n = 0
          for q in range(0, numQuadPts):
            for p in range(0, numQuadPts):
              self.quadPts[n][0] = quadpts[p]
              self.quadPts[n][1] = quadpts[q]
              self.quadWts[n]    = quadwts[p]*quadwts[q]
              
              m = 0
              for (bp,bq) in vertexOrder:
                self.basis[n][m] = basis[p][bp]*basis[q][bq]
                self.basisDeriv[n][m][0] = basisDeriv[p][bp][0]*basis[q][bq]
                self.basisDeriv[n][m][1] = basis[p][bp]*basisDeriv[q][bq][0]
                m += 1
              if not m == self.numCorners:
                raise RuntimeError('Invalid 2-D basis tabulation: '+str(m)+ \
                                 ' should be '+str(self.numCorners))
              n += 1
          if not n == self.numQuadPts:
            raise RuntimeError('Invalid 2-D quadrature: '+str(n)+ \
                               ' should be '+str(self.numQuadPts))

      elif dim == 3:
        # Set order of vertices and basis functions.
        # Corners
        vertexOrder = [(0,0,0), (0,1,0), (1,1,0), (1,0,0),
                       (0,0,1), (1,0,1), (1,1,1), (0,1,1)]
        # Edges
        #   Bottom front
        p = numpy.arange(2, numBasis, dtype=numpy.int32)
        q = numpy.zeros(numBasis-2, dtype=numpy.int32)
        r = numpy.zeros(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Bottom right
        p = numpy.ones(numBasis-2, dtype=numpy.int32)
        q = numpy.arange(2, numBasis, dtype=numpy.int32)
        r = numpy.zeros(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Bottom back
        p = numpy.arange(numBasis-1, 1, step=-1, dtype=numpy.int32)
        q = numpy.ones(numBasis-2, dtype=numpy.int32)
        r = numpy.zeros(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Bottom left
        p = numpy.zeros(numBasis-2, dtype=numpy.int32)
        q = numpy.arange(numBasis-1, 1, step=-1, dtype=numpy.int32)
        r = numpy.zeros(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Top front
        p = numpy.arange(2, numBasis, dtype=numpy.int32)
        q = numpy.zeros(numBasis-2, dtype=numpy.int32)
        r = numpy.ones(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Top right
        p = numpy.ones(numBasis-2, dtype=numpy.int32)
        q = numpy.arange(2, numBasis, dtype=numpy.int32)
        r = numpy.ones(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Top back
        p = numpy.arange(numBasis-1, 1, step=-1, dtype=numpy.int32)
        q = numpy.ones(numBasis-2, dtype=numpy.int32)
        r = numpy.ones(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Top left
        p = numpy.zeros(numBasis-2, dtype=numpy.int32)
        q = numpy.arange(numBasis-1, 1, step=-1, dtype=numpy.int32)
        r = numpy.ones(numBasis-2, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Middle left front
        p = numpy.zeros(numBasis-2, dtype=numpy.int32)
        q = numpy.zeros(numBasis-2, dtype=numpy.int32)
        r = numpy.arange(2, numBasis, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Middle right front
        p = numpy.ones(numBasis-2, dtype=numpy.int32)
        q = numpy.zeros(numBasis-2, dtype=numpy.int32)
        r = numpy.arange(2, numBasis, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Middle right back
        p = numpy.ones(numBasis-2, dtype=numpy.int32)
        q = numpy.ones(numBasis-2, dtype=numpy.int32)
        r = numpy.arange(2, numBasis, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        #   Middle left back
        p = numpy.zeros(numBasis-2, dtype=numpy.int32)
        q = numpy.ones(numBasis-2, dtype=numpy.int32)
        r = numpy.arange(2, numBasis, dtype=numpy.int32)
        vertexOrder += zip(p,q,r)
        # Face
        if numBasis > 2:
          # Left / Right
          ip = numpy.arange(0, 2, dtype=numpy.int32)
          p = numpy.tile(ip, ((numBasis-2)*(numBasis-2), 1)).transpose()
          iq = numpy.arange(2, numBasis, dtype=numpy.int32)
          q = numpy.tile(iq, (1, 2*(numBasis-2)))
          ir = numpy.arange(2, numBasis, dtype=numpy.int32)
          r = numpy.tile(ir, (2, numBasis-2)).transpose()
          vertexOrder += zip(p.ravel(),q.ravel(),r.ravel())
          # Front / Back
          ip = numpy.arange(2, numBasis, dtype=numpy.int32)
          p = numpy.tile(ip, (1, 2*(numBasis-2)))
          iq = numpy.arange(0, 2, dtype=numpy.int32)
          q = numpy.tile(iq, ((numBasis-2)*(numBasis-2), 1)).transpose()
          ir = numpy.arange(2, numBasis, dtype=numpy.int32)
          r = numpy.tile(ir, (2, numBasis-2)).transpose()
          vertexOrder += zip(p.ravel(),q.ravel(),r.ravel())
          # Bottom / Top
          ip = numpy.arange(2, numBasis, dtype=numpy.int32)
          p = numpy.tile(ip, (1, 2*(numBasis-2)))
          iq = numpy.arange(2, numBasis, dtype=numpy.int32)
          q = numpy.tile(iq, (2, numBasis-2)).transpose()
          ir = numpy.arange(0, 2, dtype=numpy.int32)
          r = numpy.tile(ir, ((numBasis-2)*(numBasis-2), 1)).transpose()
          vertexOrder += zip(p.ravel(),q.ravel(),r.ravel())          

          # Interior
          ip = numpy.arange(2, numBasis, dtype=numpy.int32)
          p = numpy.tile(ip, (1, (numBasis-2)*(numBasis-2)))
          iq = numpy.arange(2, numBasis, dtype=numpy.int32)
          q = numpy.tile(iq, ((numBasis-2), numBasis-2)).transpose()
          ir = numpy.arange(2, numBasis, dtype=numpy.int32)
          r = numpy.tile(ir, ((numBasis-2)*(numBasis-2), 1)).transpose()
          vertexOrder += zip(p.ravel(),q.ravel(),r.ravel())
        
        self.vertices = numpy.zeros((self.numCorners, dim))
        n = 0
        for (p,q,r) in vertexOrder:
          self.vertices[n][0] = vertices[p]
          self.vertices[n][1] = vertices[q]
          self.vertices[n][2] = vertices[r]
          n += 1
        if not n == self.numCorners:
          raise RuntimeError('Invalid 3-D vertex ordering: '+str(n)+ \
                             ' should be '+str(self.numCorners))

        self.quadPts    = numpy.zeros((numQuadPts*numQuadPts*numQuadPts, dim))
        self.quadWts    = numpy.zeros((numQuadPts*numQuadPts*numQuadPts,))
        self.basis      = numpy.zeros((numQuadPts*numQuadPts*numQuadPts,
                                       numBasis*numBasis*numBasis))
        self.basisDeriv = numpy.zeros((numQuadPts*numQuadPts*numQuadPts,
                                       numBasis*numBasis*numBasis,
                                       dim))

        # Order of quadrature points doesn't matter
        # Order of basis functions should match vertices for isoparametric
        n = 0
        for r in range(0, numQuadPts):
          for q in range(0, numQuadPts):
            for p in range(0, numQuadPts):
              self.quadPts[n][0] = quadpts[p]
              self.quadPts[n][1] = quadpts[q]
              self.quadPts[n][2] = quadpts[r]
              self.quadWts[n]    = quadwts[p]*quadwts[q]*quadwts[r]
            
              m = 0
              for (bp,bq,br) in vertexOrder:
                self.basis[n][m] = basis[p][bp]*basis[q][bq]*basis[r][br]
                self.basisDeriv[n][m][0] = basisDeriv[p][bp][0]*basis[q][bq]*basis[r][br]
                self.basisDeriv[n][m][1] = basis[p][bp]*basisDeriv[q][bq][0]*basis[r][br]
                self.basisDeriv[n][m][2] = basis[p][bp]*basis[q][bq]*basisDeriv[r][br][0]
                m += 1

              if not m == self.numCorners:
                raise RuntimeError('Invalid 3-D basis tabulation: '+str(m)+ \
                                   ' should be '+str(self.numCorners))
              n += 1
        if not n == self.numQuadPts:
          raise RuntimeError('Invalid 3-D quadrature: '+str(n)+ \
                             ' should be '+str(self.numQuadPts))

      self.vertices = numpy.reshape(self.vertices, (self.numCorners, dim))
      self.quadPts = numpy.reshape(self.quadPts, (self.numQuadPts, dim))
      self.quadWts = numpy.reshape(self.quadWts, (self.numQuadPts))
      self.basis = numpy.reshape(self.basis, (self.numQuadPts, self.numCorners))
      self.basisDeriv = numpy.reshape(self.basisDeriv, (self.numQuadPts, self.numCorners, dim))
    else:
      # Need 0-D quadrature for boundary conditions of 1-D meshes
      self.cellDim = 0
      self.numCorners = 1
      self.numQuadPts = 1
      self.basis = numpy.array([1.0], dtype=numpy.float64)
      self.basisDeriv = numpy.array([1.0], dtype=numpy.float64)
      self.quadPts = numpy.array([0.0], dtype=numpy.float64)
      self.quadWts = numpy.array([1.0], dtype=numpy.float64)

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
        self.order = self.degree+1
    except ValueError, err:
      aliases = ", ".join(self.aliases)
      raise ValueError("Error while configuring quadrature "
                       "(%s):\n%s" % (aliases, err.message))
    return


  def _setupGeometry(self, spaceDim):
    """
    Setup reference cell geometry object.
    """
    import CellGeometry
    self.geometry = None
    if 3 == self.cellDim:
      if 3 == spaceDim:
        self.geometry = CellGeometry.GeometryHex3D()
    elif 2 == self.cellDim:
      if 2 == spaceDim:
        self.geometry = CellGeometry.GeometryQuad2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryQuad3D()
    elif 1 == self.cellDim:
      if 2 == spaceDim:
        self.geometry = CellGeometry.GeometryLine2D()
      elif 3 == spaceDim:
        self.geometry = CellGeometry.GeometryLine3D()
    if None == self.geometry:
      raise ValueError("Could not set shape of %dD cell for "
                       "Quadrature object '%s' in spatial " \
                       "dimension '%d'." % (self.cellDim, self.name, spaceDim))
    return
  

  def _setupQuadrature(self):
    """
    Setup quadrature rule for reference cell.
    """
    from FIAT.reference_element import default_simplex
    from FIAT.quadrature import make_quadrature
    from FIATQuadrature import CollocatedQuadratureRule
      
    if not self.collocateQuad:
      q = make_quadrature(default_simplex(1), self.order)
    else:
      q = CollocatedQuadratureRule(default_simplex(1), self.order)

    return q


  def _setupElement(self):
    """
    Setup the finite element for reference cell.
    """
    from FIAT.lagrange import Lagrange
    from FIAT.reference_element import default_simplex
    return Lagrange(default_simplex(1), self.degree)


  def _setupVertices(self, element):
    """
    Setup evaluation functional points for reference cell.
    """
    return [n.get_point_dict().keys()[0] for n in element.dual.get_nodes()]


# FACTORIES ////////////////////////////////////////////////////////////

def reference_cell():
  """
  Factory associated with FIATLagrange.
  """
  return FIATLagrange()


# End of file 

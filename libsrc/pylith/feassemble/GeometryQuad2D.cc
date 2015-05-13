// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "GeometryQuad2D.hh" // implementation of class methods

#include "GeometryLine2D.hh" // USES GeometryLine2D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES scalar_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryQuad2D::GeometryQuad2D(void) :
  CellGeometry(QUADRILATERAL, 2)
{ // constructor
  const PylithScalar vertices[4*2] = {
    -1.0,  -1.0,
    +1.0,  -1.0,
    +1.0,  +1.0,
    -1.0,  +1.0,
  };
  _setVertices(vertices, 4, 2);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryQuad2D::~GeometryQuad2D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryQuad2D::clone(void) const
{ // clone
  return new GeometryQuad2D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryQuad2D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryLine2D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryQuad2D::ptsRefToGlobal(PylithScalar* ptsGlobal,
						   const PylithScalar* ptsRef,
						   const PylithScalar* vertices,
						   const int dim,
						   const int npts) const
{ // ptsRefToGlobal
  assert(ptsGlobal);
  assert(ptsRef);
  assert(vertices);
  assert(2 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];

  const PylithScalar x1 = vertices[2];
  const PylithScalar y1 = vertices[3];

  const PylithScalar x2 = vertices[4];
  const PylithScalar y2 = vertices[5];

  const PylithScalar x3 = vertices[6];
  const PylithScalar y3 = vertices[7];

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;

  const PylithScalar f_3 = x3 - x0;
  const PylithScalar g_3 = y3 - y0;

  const PylithScalar f_01 = x2 - x1 - x3 + x0;
  const PylithScalar g_01 = y2 - y1 - y3 + y0;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p1 = 0.5 * (1.0 + ptsRef[iR++]);
    ptsGlobal[iG++] = x0 + f_1 * p0 + f_3 * p1 + f_01 * p0 * p1;
    ptsGlobal[iG++] = y0 + g_1 * p0 + g_3 * p1 + g_01 * p0 * p1;
  } // for

  PetscLogFlops(10 + npts*18);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad2D::jacobian(scalar_array* jacobian,
					     PylithScalar* det,
					     const PylithScalar* vertices,
					     const int numVertices,
					     const int spaceDim,
					     const PylithScalar* location,
					     const int cellDim) const
{ // jacobian
  assert(jacobian);
  assert(det);
  assert(vertices);
  assert(location);

  assert(this->numCorners() == numVertices || // linear
	 this->numCorners()+1 == numVertices); // quadratic
  assert(this->spaceDim() == spaceDim);
  assert(this->cellDim() == cellDim);
  assert(size_t(spaceDim*cellDim) == jacobian->size());
  
  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];

  const PylithScalar x1 = vertices[2];
  const PylithScalar y1 = vertices[3];

  const PylithScalar x2 = vertices[4];
  const PylithScalar y2 = vertices[5];

  const PylithScalar x3 = vertices[6];
  const PylithScalar y3 = vertices[7];

  const PylithScalar x = 0.5 * (location[0] + 1.0);
  const PylithScalar y = 0.5 * (location[1] + 1.0);
  assert(0 <= x && x <= 1.0);
  assert(0 <= y && y <= 1.0);

  const PylithScalar f_xy = x2 - x1 - x3 + x0;
  const PylithScalar g_xy = y2 - y1 - y3 + y0;

  (*jacobian)[0] = (x1 - x0 + f_xy*y) / 2.0;
  (*jacobian)[1] = (x3 - x0 + f_xy*x) / 2.0;
  (*jacobian)[2] = (y1 - y0 + g_xy*y) / 2.0;
  (*jacobian)[3] = (y3 - y0 + g_xy*x) / 2.0;

  *det = 
    (*jacobian)[0]*(*jacobian)[3] - 
    (*jacobian)[1]*(*jacobian)[2];

  PetscLogFlops(31);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad2D::jacobian(PylithScalar* jacobian,
					     PylithScalar* det,
					     const PylithScalar* vertices,
					     const PylithScalar* ptsRef,
					     const int dim,
					     const int npts) const
{ // jacobian
  assert(jacobian);
  assert(det);
  assert(vertices);
  assert(ptsRef);
  assert(2 == dim);
  assert(spaceDim() == dim);
    
  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];

  const PylithScalar x1 = vertices[2];
  const PylithScalar y1 = vertices[3];

  const PylithScalar x2 = vertices[4];
  const PylithScalar y2 = vertices[5];

  const PylithScalar x3 = vertices[6];
  const PylithScalar y3 = vertices[7];

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;

  const PylithScalar f_3 = x3 - x0;
  const PylithScalar g_3 = y3 - y0;

  const PylithScalar f_01 = x2 - x1 - x3 + x0;
  const PylithScalar g_01 = y2 - y1 - y3 + y0;

  for (int i=0, iR=0, iJ=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p1 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    const PylithScalar j00 = (f_1 + f_01 * p1) / 2.0; 
    const PylithScalar j01 = (f_3 + f_01 * p0) / 2.0; 
    const PylithScalar j10 = (g_1 + g_01 * p1) / 2.0;
    const PylithScalar j11 = (g_3 + g_01 * p0) / 2.0; 

    jacobian[iJ++] = j00;
    jacobian[iJ++] = j01;
    jacobian[iJ++] = j10;
    jacobian[iJ++] = j11;
    det[i] = j00*j11 - j01*j10;
  } // for

  PetscLogFlops(10 + npts*19);
} // jacobian


// ----------------------------------------------------------------------
// Compute minimum width across cell.
PylithScalar
pylith::feassemble::GeometryQuad2D::minCellWidth(const PylithScalar* coordinatesCell,
						 const int numVertices,
						 const int spaceDim) const
{ // minCellWidth
  const int numCorners = 4;
  const int dim = 2;

  assert(numCorners == numVertices);
  assert(dim == spaceDim);

  const int numEdges = 4;
  const int edges[numEdges][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},
  };

  PylithScalar minWidth = PYLITH_MAXSCALAR;
  for (int iedge=0; iedge < numEdges; ++iedge) {
    const int iA = edges[iedge][0];
    const int iB = edges[iedge][1];
    const PylithScalar xA = coordinatesCell[dim*iA  ];
    const PylithScalar yA = coordinatesCell[dim*iA+1];
    const PylithScalar xB = coordinatesCell[dim*iB  ];
    const PylithScalar yB = coordinatesCell[dim*iB+1];
    
    const PylithScalar edgeLen = sqrt(pow(xB-xA,2) + pow(yB-yA,2));
    if (edgeLen < minWidth) {
      minWidth = edgeLen;
    } // if
  } // for

  PetscLogFlops(numEdges*6);

  return minWidth;
} // minCellWidth


// End of file

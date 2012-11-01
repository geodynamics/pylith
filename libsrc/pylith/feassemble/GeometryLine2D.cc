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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "GeometryLine2D.hh" // implementation of class methods

#include "GeometryPoint2D.hh" // USES GeometryPoint

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES scalar_array

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryLine2D::GeometryLine2D(void) :
  CellGeometry(LINE, 2)
{ // constructor
  const PylithScalar vertices[] = {
    -1.0,
    +1.0,
  };
  _setVertices(vertices, 2, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryLine2D::~GeometryLine2D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine2D::clone(void) const
{ // clone
  return new GeometryLine2D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine2D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryPoint2D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryLine2D::ptsRefToGlobal(PylithScalar* ptsGlobal,
						   const PylithScalar* ptsRef,
						   const PylithScalar* vertices,
						   const int dim,
						   const int npts) const
{ // ptsRefToGlobal
  assert(0 != ptsGlobal);
  assert(0 != ptsRef);
  assert(0 != vertices);
  assert(2 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];

  const PylithScalar x1 = vertices[2];
  const PylithScalar y1 = vertices[3];

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;

  for (int i=0, iG=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[i]);
    ptsGlobal[iG++] = x0 + f_1 * p0;
    ptsGlobal[iG++] = y0 + g_1 * p0;
  } // for

  PetscLogFlops(2 + npts*6);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine2D::jacobian(scalar_array* jacobian,
					     PylithScalar* det,
					     const scalar_array& vertices,
					     const scalar_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert( (numCorners()*spaceDim() == vertices.size()) || // linear
	  ((numCorners()+1)*spaceDim() == vertices.size()) ); // quadratic
  assert(spaceDim()*cellDim() == jacobian->size());

  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];

  const PylithScalar x1 = vertices[2];
  const PylithScalar y1 = vertices[3];

  (*jacobian)[0] = (x1 - x0)/2.0;
  (*jacobian)[1] = (y1 - y0)/2.0;
  *det = sqrt(pow((*jacobian)[0], 2) +
	      pow((*jacobian)[1], 2));

  PetscLogFlops(8);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine2D::jacobian(PylithScalar* jacobian,
					     PylithScalar* det,
					     const PylithScalar* vertices,
					     const PylithScalar* location,
					     const int dim,
					     const int npts) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != location);
  assert(2 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];

  const PylithScalar x1 = vertices[2];
  const PylithScalar y1 = vertices[3];

  const PylithScalar j1 = (x1 - x0) / 2.0;
  const PylithScalar j2 = (y1 - y0) / 2.0;
  const PylithScalar jdet = sqrt(j1*j1 + j2*j2);

  for (int i=0, iJ=0; i < npts; ++i) {
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    det[i] = jdet;
  } // for

  PetscLogFlops(8);
} // jacobian


// ----------------------------------------------------------------------
// Compute minimum width across cell.
PylithScalar
pylith::feassemble::GeometryLine2D::minCellWidth(const scalar_array& coordinatesCell) const
{ // minCellWidth
  const int numCorners = 2;
  const int spaceDim = 2;
  assert(numCorners*spaceDim == coordinatesCell.size());

  const PylithScalar xA = coordinatesCell[0];
  const PylithScalar yA = coordinatesCell[1];
  const PylithScalar xB = coordinatesCell[2];
  const PylithScalar yB = coordinatesCell[3];
    
  const PylithScalar minWidth = sqrt(pow(xB-xA,2) + pow(yB-yA,2));

  PetscLogFlops(6);

  return minWidth;
} // minCellWidth


// End of file

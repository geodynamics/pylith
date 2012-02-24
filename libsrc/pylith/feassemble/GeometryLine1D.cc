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

#include "GeometryLine1D.hh" // implementation of class methods

#include "GeometryPoint1D.hh" // USES GeometryPoint

#include "pylith/utils/array.hh" // USES scalar_array

#include "petsc.h" // USES PetscLogFlops

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryLine1D::GeometryLine1D(void) :
  CellGeometry(LINE, 1)
{ // constructor
  const PylithScalar vertices[] = {
    -1.0,
    +1.0,
  };
  _setVertices(vertices, 2, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryLine1D::~GeometryLine1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine1D::clone(void) const
{ // clone
  return new GeometryLine1D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine1D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryPoint1D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryLine1D::ptsRefToGlobal(PylithScalar* ptsGlobal,
						   const PylithScalar* ptsRef,
						   const PylithScalar* vertices,
						   const int dim,
						   const int npts) const
{ // ptsRefToGlobal
  assert(0 != ptsGlobal);
  assert(0 != ptsRef);
  assert(0 != vertices);
  assert(1 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar x1 = vertices[1];

  for (int i=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5*(1.0+ptsRef[i]);
    ptsGlobal[i] = x0 + (x1-x0) * p0;
  } // for

  PetscLogFlops(5*npts);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine1D::jacobian(scalar_array* jacobian,
					     PylithScalar* det,
					     const scalar_array& vertices,
					     const scalar_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert( (numCorners()*spaceDim() == vertices.size()) || // linear edge
	  ((numCorners()+1)*spaceDim() == vertices.size()) ); // quadratic edge
  assert(spaceDim()*cellDim() == jacobian->size());

  const PylithScalar x0 = vertices[0];
  const PylithScalar x1 = vertices[1];

  (*jacobian)[0] = (x1 - x0)/2.0;
  *det = (*jacobian)[0];
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine1D::jacobian(PylithScalar* jacobian,
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
  assert(1 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar x1 = vertices[1];

  const PylithScalar j = (x1 - x0)/2.0;
  for (int i=0; i < npts; ++i) {
    jacobian[i] = j;
    det[i] = j;
  } // for

  PetscLogFlops(2);
} // jacobian


// End of file

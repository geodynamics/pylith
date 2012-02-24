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

#include "GeometryLine3D.hh" // implementation of class methods

#include "GeometryPoint3D.hh" // USES GeometryPoint3D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES double_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryLine3D::GeometryLine3D(void) :
  CellGeometry(LINE, 3)
{ // constructor
  const double vertices[] = {
    -1.0,
    +1.0,
  };
  _setVertices(vertices, 2, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryLine3D::~GeometryLine3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine3D::clone(void) const
{ // clone
  return new GeometryLine3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryPoint3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryLine3D::ptsRefToGlobal(double* ptsGlobal,
						   const double* ptsRef,
						   const double* vertices,
						   const int dim,
						   const int npts) const
{ // ptsRefToGlobal
  assert(0 != ptsGlobal);
  assert(0 != ptsRef);
  assert(0 != vertices);
  assert(3 == dim);
  assert(spaceDim() == dim);

  const double x0 = vertices[0];
  const double y0 = vertices[1];
  const double z0 = vertices[2];

  const double x1 = vertices[3];
  const double y1 = vertices[4];
  const double z1 = vertices[5];

  const double f_1 = x1 - x0;
  const double g_1 = y1 - y0;
  const double h_1 = z1 - z0;

  for (int i=0, iG=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[i]);
    ptsGlobal[iG++] = x0 + f_1 * p0;
    ptsGlobal[iG++] = y0 + g_1 * p0;
    ptsGlobal[iG++] = z0 + h_1 * p0;
  } // for

  PetscLogFlops(3 + npts*8);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine3D::jacobian(double_array* jacobian,
					   double* det,
					   const double_array& vertices,
					   const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert( (numCorners()*spaceDim() == vertices.size()) || // linear edge
	  ((numCorners()+1)*spaceDim() == vertices.size()) ); // quadratic edge
  assert(spaceDim()*cellDim() == jacobian->size());

  const double x0 = vertices[0];
  const double y0 = vertices[1];
  const double z0 = vertices[2];

  const double x1 = vertices[3];
  const double y1 = vertices[4];
  const double z1 = vertices[5];

  (*jacobian)[0] = (x1 - x0)/2.0;
  (*jacobian)[1] = (y1 - y0)/2.0;
  (*jacobian)[2] = (z1 - z0)/2.0;
  *det = sqrt(pow((*jacobian)[0], 2) +
	      pow((*jacobian)[1], 2) +
	      pow((*jacobian)[2], 2));
  PetscLogFlops(12);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine3D::jacobian(double* jacobian,
					     double* det,
					     const double* vertices,
					     const double* location,
					     const int dim,
					     const int npts) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != location);
  assert(3 == dim);
  assert(spaceDim() == dim);

  const double x0 = vertices[0];
  const double y0 = vertices[1];
  const double z0 = vertices[2];

  const double x1 = vertices[3];
  const double y1 = vertices[4];
  const double z1 = vertices[5];

  const double j1 = (x1 - x0) / 2.0;
  const double j2 = (y1 - y0) / 2.0;
  const double j3 = (z1 - z0) / 2.0;
  const double jdet = sqrt(j1*j1 + j2*j2 + j3*j3);

  for (int i=0, iJ=0; i < npts; ++i) {
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    jacobian[iJ++] = j3;
    det[i] = jdet;
  } // for

  PetscLogFlops(12);
} // jacobian


// End of file

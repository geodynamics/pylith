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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "GeometryQuad3D.hh" // implementation of class methods

#include "GeometryLine3D.hh" // USES GeometryLine3D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES double_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryQuad3D::GeometryQuad3D(void) :
  CellGeometry(QUADRILATERAL, 3)
{ // constructor
  const double vertices[] = {
    -1.0,  -1.0,
    +1.0,  -1.0,
    +1.0,  +1.0,
    -1.0,  +1.0,
  };
  _setVertices(vertices, 4, 2);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryQuad3D::~GeometryQuad3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryQuad3D::clone(void) const
{ // clone
  return new GeometryQuad3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryQuad3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryLine3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryQuad3D::ptsRefToGlobal(double* ptsGlobal,
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

  const double x2 = vertices[6];
  const double y2 = vertices[7];
  const double z2 = vertices[8];

  const double x3 = vertices[9];
  const double y3 = vertices[10];
  const double z3 = vertices[11];

  const double f_1 = x1 - x0;
  const double g_1 = y1 - y0;
  const double h_1 = z1 - z0;

  const double f_3 = x3 - x0;
  const double g_3 = y3 - y0;
  const double h_3 = z3 - z0;

  const double f_01 = x2 - x1 - x3 + x0;
  const double g_01 = y2 - y1 - y3 + y0;
  const double h_01 = z2 - z1 - z3 + z0;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p1 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    ptsGlobal[iG++] = x0 + f_1 * p0 + f_3 * p1 + f_01 * p0 * p1;
    ptsGlobal[iG++] = y0 + g_1 * p0 + g_3 * p1 + g_01 * p0 * p1;
    ptsGlobal[iG++] = z0 + h_1 * p0 + h_3 * p1 + h_01 * p0 * p1;
  } // for

  PetscLogFlops(15 + npts*25);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad3D::jacobian(double_array* jacobian,
					  double* det,
					  const double_array& vertices,
					  const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert( (numCorners()*spaceDim() == vertices.size()) || // linear quad
	  ((numCorners()+5)*spaceDim() == vertices.size()) ); // quadratic quad
  assert(cellDim() == location.size());
  assert(spaceDim()*cellDim() == jacobian->size());
  
  const double x0 = vertices[0];
  const double y0 = vertices[1];
  const double z0 = vertices[2];

  const double x1 = vertices[3];
  const double y1 = vertices[4];
  const double z1 = vertices[5];

  const double x2 = vertices[6];
  const double y2 = vertices[7];
  const double z2 = vertices[8];

  const double x3 = vertices[9];
  const double y3 = vertices[10];
  const double z3 = vertices[11];

  const double x = 0.5 * (location[0] + 1.0);
  const double y = 0.5 * (location[1] + 1.0);
  assert(0 <= x && x <= 1.0);
  assert(0 <= y && y <= 1.0);

  const double f_xy = x2 - x1 - x3 + x0;
  const double g_xy = y2 - y1 - y3 + y0;
  const double h_xy = z2 - z1 - z3 + z0;

  (*jacobian)[0] = (x1 - x0 + f_xy*y) / 2.0;
  (*jacobian)[1] = (x3 - x0 + f_xy*x) / 2.0;

  (*jacobian)[2] = (y1 - y0 + g_xy*y) / 2.0;
  (*jacobian)[3] = (y3 - y0 + g_xy*x) / 2.0;

  (*jacobian)[4] = (z1 - z0 + h_xy*y) / 2.0;
  (*jacobian)[5] = (z3 - z0 + h_xy*x) / 2.0;

  const double jj00 = 
    (*jacobian)[0]*(*jacobian)[0] +
    (*jacobian)[2]*(*jacobian)[2] +
    (*jacobian)[4]*(*jacobian)[4];
  const double jj10 =
    (*jacobian)[0]*(*jacobian)[1] +
    (*jacobian)[2]*(*jacobian)[3] +
    (*jacobian)[4]*(*jacobian)[5];
  const double jj01 = jj10;
  const double jj11 = 
    (*jacobian)[1]*(*jacobian)[1] +
    (*jacobian)[3]*(*jacobian)[3] +
    (*jacobian)[5]*(*jacobian)[5];
  *det = sqrt(jj00*jj11 - jj01*jj10);

  PetscLogFlops(50);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad3D::jacobian(double* jacobian,
					     double* det,
					     const double* vertices,
					     const double* ptsRef,
					     const int dim,
					     const int npts) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != ptsRef);
  assert(3 == dim);
  assert(spaceDim() == dim);
  
  const double x0 = vertices[0];
  const double y0 = vertices[1];
  const double z0 = vertices[2];

  const double x1 = vertices[3];
  const double y1 = vertices[4];
  const double z1 = vertices[5];

  const double x2 = vertices[6];
  const double y2 = vertices[7];
  const double z2 = vertices[8];

  const double x3 = vertices[9];
  const double y3 = vertices[10];
  const double z3 = vertices[11];

  const double f_1 = (x1 - x0) / 2.0;
  const double g_1 = (y1 - y0) / 2.0;
  const double h_1 = (z1 - z0) / 2.0;

  const double f_3 = (x3 - x0) / 2.0;
  const double g_3 = (y3 - y0) / 2.0;
  const double h_3 = (z3 - z0) / 2.0;

  const double f_01 = (x2 - x1 - x3 + x0) / 2.0;
  const double g_01 = (y2 - y1 - y3 + y0) / 2.0;
  const double h_01 = (z2 - z1 - z3 + z0) / 2.0;

  for (int i=0, iR=0, iJ=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p1 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    const double j0 = f_1 + f_01 * p1; 
    const double j1 = f_3 + f_01 * p0; 
    const double j2 = g_1 + g_01 * p1;
    const double j3 = g_3 + g_01 * p0; 
    const double j4 = h_1 + h_01 * p1;
    const double j5 = h_3 + h_01 * p0;
    jacobian[iJ++] = j0;
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    jacobian[iJ++] = j3;
    jacobian[iJ++] = j4;
    jacobian[iJ++] = j5;

    const double jj00 = j0*j0 + j2*j2 + j4*j4;
    const double jj10 = j0*j1 + j2*j3 + j4*j5;
    const double jj01 = jj10;
    const double jj11 = j1*j1 + j3*j3 + j5*j5;
    det[i] = sqrt(jj00*jj11 - jj01*jj10);
  } // for

  PetscLogFlops(28 + npts*32);
} // jacobian


// End of file

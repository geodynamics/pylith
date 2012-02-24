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

#include "GeometryTet3D.hh" // implementation of class methods

#include "GeometryTri3D.hh" // USES GeometryTri3D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES scalar_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryTet3D::GeometryTet3D(void) :
  CellGeometry(TETRAHEDRON, 3)
{ // constructor
  const PylithScalar vertices[] = {
    -1.0,  -1.0,  -1.0,
    +1.0,  -1.0,  -1.0,
    -1.0,  +1.0,  -1.0,
    -1.0,  -1.0,  +1.0,
  };
  _setVertices(vertices, 4, 3);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryTet3D::~GeometryTet3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTet3D::clone(void) const
{ // clone
  return new GeometryTet3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTet3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryTri3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryTet3D::ptsRefToGlobal(PylithScalar* ptsGlobal,
						  const PylithScalar* ptsRef,
						  const PylithScalar* vertices,
						  const int dim,
						  const int npts) const
{ // ptsRefToGlobal
  assert(0 != ptsGlobal);
  assert(0 != ptsRef);
  assert(0 != vertices);
  assert(3 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];
  const PylithScalar z0 = vertices[2];

  const PylithScalar x1 = vertices[3];
  const PylithScalar y1 = vertices[4];
  const PylithScalar z1 = vertices[5];

  const PylithScalar x2 = vertices[6];
  const PylithScalar y2 = vertices[7];
  const PylithScalar z2 = vertices[8];

  const PylithScalar x3 = vertices[9];
  const PylithScalar y3 = vertices[10];
  const PylithScalar z3 = vertices[11];

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;
  const PylithScalar h_1 = z1 - z0;

  const PylithScalar f_2 = x2 - x0;
  const PylithScalar g_2 = y2 - y0;
  const PylithScalar h_2 = z2 - z0;

  const PylithScalar f_3 = x3 - x0;
  const PylithScalar g_3 = y3 - y0;
  const PylithScalar h_3 = z3 - z0;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p1 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p2 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    assert(0 <= p2 && p2 <= 1.0);

    ptsGlobal[iG++] = x0 + f_1 * p0 + f_2 * p1 + f_3 * p2;
    ptsGlobal[iG++] = y0 + g_1 * p0 + g_2 * p1 + g_3 * p2;;
    ptsGlobal[iG++] = z0 + h_1 * p0 + h_2 * p1 + h_3 * p2;
  } // for

  PetscLogFlops(9 + npts*24);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTet3D::jacobian(scalar_array* jacobian,
					    PylithScalar* det,
					    const scalar_array& vertices,
					    const scalar_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert(numCorners()*spaceDim() == vertices.size());
  assert(spaceDim()*cellDim() == jacobian->size());
  
  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];
  const PylithScalar z0 = vertices[2];

  const PylithScalar x1 = vertices[3];
  const PylithScalar y1 = vertices[4];
  const PylithScalar z1 = vertices[5];

  const PylithScalar x2 = vertices[6];
  const PylithScalar y2 = vertices[7];
  const PylithScalar z2 = vertices[8];

  const PylithScalar x3 = vertices[9];
  const PylithScalar y3 = vertices[10];
  const PylithScalar z3 = vertices[11];

  (*jacobian)[0] = (x1 - x0) / 2.0;
  (*jacobian)[1] = (x2 - x0) / 2.0;
  (*jacobian)[2] = (x3 - x0) / 2.0;
  (*jacobian)[3] = (y1 - y0) / 2.0;
  (*jacobian)[4] = (y2 - y0) / 2.0;
  (*jacobian)[5] = (y3 - y0) / 2.0;
  (*jacobian)[6] = (z1 - z0) / 2.0;
  (*jacobian)[7] = (z2 - z0) / 2.0;
  (*jacobian)[8] = (z3 - z0) / 2.0;

  *det = 
    (*jacobian)[0]*((*jacobian)[4]*(*jacobian)[8] -
		    (*jacobian)[5]*(*jacobian)[7]) -
    (*jacobian)[1]*((*jacobian)[3]*(*jacobian)[8] -
		    (*jacobian)[5]*(*jacobian)[6]) +
    (*jacobian)[2]*((*jacobian)[3]*(*jacobian)[7] -
		    (*jacobian)[4]*(*jacobian)[6]);

  PetscLogFlops(32);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTet3D::jacobian(PylithScalar* jacobian,
					    PylithScalar* det,
					    const PylithScalar* vertices,
					    const PylithScalar* ptsRef,
					    const int dim,
					    const int npts) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != ptsRef);
  assert(3 == dim);
  assert(spaceDim() == dim);
  
  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];
  const PylithScalar z0 = vertices[2];

  const PylithScalar x1 = vertices[3];
  const PylithScalar y1 = vertices[4];
  const PylithScalar z1 = vertices[5];

  const PylithScalar x2 = vertices[6];
  const PylithScalar y2 = vertices[7];
  const PylithScalar z2 = vertices[8];

  const PylithScalar x3 = vertices[9];
  const PylithScalar y3 = vertices[10];
  const PylithScalar z3 = vertices[11];

  const PylithScalar j0 = (x1 - x0) / 2.0;
  const PylithScalar j1 = (x2 - x0) / 2.0;
  const PylithScalar j2 = (x3 - x0) / 2.0;

  const PylithScalar j3 = (y1 - y0) / 2.0;
  const PylithScalar j4 = (y2 - y0) / 2.0;
  const PylithScalar j5 = (y3 - y0) / 2.0;

  const PylithScalar j6 = (z1 - z0) / 2.0;
  const PylithScalar j7 = (z2 - z0) / 2.0;
  const PylithScalar j8 = (z3 - z0) / 2.0;

  const PylithScalar jdet = 
    j0*(j4*j8 - j5*j7) -
    j1*(j3*j8 - j5*j6) +
    j2*(j3*j7 - j4*j6);

  for (int i=0, iJ=0; i < npts; ++i) {
    jacobian[iJ++] = j0;
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    jacobian[iJ++] = j3;
    jacobian[iJ++] = j4;
    jacobian[iJ++] = j5;
    jacobian[iJ++] = j6;
    jacobian[iJ++] = j7;
    jacobian[iJ++] = j8;
    det[i] = jdet;
  } // for

  PetscLogFlops(32);
} // jacobian


// End of file

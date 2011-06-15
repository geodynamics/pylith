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

#include "GeometryHex3D.hh" // implementation of class methods

#include "GeometryQuad3D.hh" // USES GeometryQuad3D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES double_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryHex3D::GeometryHex3D(void) :
  CellGeometry(HEXAHEDRON, 3)
{ // constructor
  const double vertices[] = {
    -1.0,  -1.0,  -1.0,
    +1.0,  -1.0,  -1.0,
    +1.0,  +1.0,  -1.0,
    -1.0,  +1.0,  -1.0,
    -1.0,  -1.0,  +1.0,
    +1.0,  -1.0,  +1.0,
    +1.0,  +1.0,  +1.0,
    -1.0,  +1.0,  +1.0,
  };
  _setVertices(vertices, 8, 3);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryHex3D::~GeometryHex3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryHex3D::clone(void) const
{ // clone
  return new GeometryHex3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryHex3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryQuad3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryHex3D::ptsRefToGlobal(double* ptsGlobal,
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

  const double x4 = vertices[12];
  const double y4 = vertices[13];
  const double z4 = vertices[14];

  const double x5 = vertices[15];
  const double y5 = vertices[16];
  const double z5 = vertices[17];

  const double x6 = vertices[18];
  const double y6 = vertices[19];
  const double z6 = vertices[20];

  const double x7 = vertices[21];
  const double y7 = vertices[22];
  const double z7 = vertices[23];

  const double f_1 = x1 - x0;
  const double g_1 = y1 - y0;
  const double h_1 = z1 - z0;

  const double f_3 = x3 - x0;
  const double g_3 = y3 - y0;
  const double h_3 = z3 - z0;

  const double f_4 = x4 - x0;
  const double g_4 = y4 - y0;
  const double h_4 = z4 - z0;

  const double f_01 = x2 - x1 - x3 + x0;
  const double g_01 = y2 - y1 - y3 + y0;
  const double h_01 = z2 - z1 - z3 + z0;

  const double f_12 = x7 - x3 - x4 + x0;
  const double g_12 = y7 - y3 - y4 + y0;
  const double h_12 = z7 - z3 - z4 + z0;

  const double f_02 = x5 - x1 - x4 + x0;
  const double g_02 = y5 - y1 - y4 + y0;
  const double h_02 = z5 - z1 - z4 + z0;

  const double f_012 = x6 - x0 + x1 - x2 + x3 + x4 - x5 - x7;
  const double g_012 = y6 - y0 + y1 - y2 + y3 + y4 - y5 - y7;
  const double h_012 = z6 - z0 + z1 - z2 + z3 + z4 - z5 - z7;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p1 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p2 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    assert(0 <= p2 && p2 <= 1.0);
    ptsGlobal[iG++] = x0 + f_1*p0 + f_3*p1 + f_4*p2
      + f_01*p0*p1 + f_12*p1*p2 + f_02*p0*p2 + f_012*p0*p1*p2;
    ptsGlobal[iG++] = y0 + g_1*p0 + g_3*p1 + g_4*p2 + g_01*p0*p1
      + g_01*p0*p1 + g_12*p1*p2 + g_02*p0*p2 + g_012*p0*p1*p2;
    ptsGlobal[iG++] = z0 + h_1*p0 + h_3*p1 + h_4*p2 + h_01*p0*p1
      + h_01*p0*p1 + h_12*p1*p2 + h_02*p0*p2 + h_012*p0*p1*p2;
  } // for

  PetscLogFlops(57 + npts*57);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryHex3D::jacobian(double_array* jacobian,
					    double* det,
					    const double_array& vertices,
					    const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert( (numCorners()*spaceDim() == vertices.size()) || // linear hex
	  ((numCorners()+19)*spaceDim() == vertices.size()) ); // quadratic hex
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

  const double x4 = vertices[12];
  const double y4 = vertices[13];
  const double z4 = vertices[14];

  const double x5 = vertices[15];
  const double y5 = vertices[16];
  const double z5 = vertices[17];

  const double x6 = vertices[18];
  const double y6 = vertices[19];
  const double z6 = vertices[20];

  const double x7 = vertices[21];
  const double y7 = vertices[22];
  const double z7 = vertices[23];

  const double x = 0.5 * (location[0] + 1.0);
  const double y = 0.5 * (location[1] + 1.0);
  const double z = 0.5 * (location[2] + 1.0);
  assert(-1.0 <= x && x <= 1.0);
  assert(-1.0 <= y && y <= 1.0);
  assert(-1.0 <= z && z <= 1.0);

  const double f_xy = x2 - x1 - x3 + x0;
  const double g_xy = y2 - y1 - y3 + y0;
  const double h_xy = z2 - z1 - z3 + z0;

  const double f_yz = x7 - x3 - x4 + x0;
  const double g_yz = y7 - y3 - y4 + y0;
  const double h_yz = z7 - z3 - z4 + z0;

  const double f_xz = x5 - x1 - x4 + x0;
  const double g_xz = y5 - y1 - y4 + y0;
  const double h_xz = z5 - z1 - z4 + z0;

  const double f_xyz = x6 - x0 + x1 - x2 + x3 + x4 - x5 - x7;
  const double g_xyz = y6 - y0 + y1 - y2 + y3 + y4 - y5 - y7;
  const double h_xyz = z6 - z0 + z1 - z2 + z3 + z4 - z5 - z7;

  (*jacobian)[0] = (x1 - x0 + f_xy*y + f_xz*z + f_xyz*y*z) / 2.0;
  (*jacobian)[1] = (x3 - x0 + f_xy*x + f_yz*z + f_xyz*x*z) / 2.0;
  (*jacobian)[2] = (x4 - x0 + f_yz*y + f_xz*x + f_xyz*x*y) / 2.0;

  (*jacobian)[3] = (y1 - y0 + g_xy*y + g_xz*z + g_xyz*y*z) / 2.0;
  (*jacobian)[4] = (y3 - y0 + g_xy*x + g_yz*z + g_xyz*x*z) / 2.0;
  (*jacobian)[5] = (y4 - y0 + g_yz*y + g_xz*x + g_xyz*x*y) / 2.0;

  (*jacobian)[6] = (z1 - z0 + h_xy*y + h_xz*z + h_xyz*y*z) / 2.0;
  (*jacobian)[7] = (z3 - z0 + h_xy*x + h_yz*z + h_xyz*x*z) / 2.0;
  (*jacobian)[8] = (z4 - z0 + h_yz*y + h_xz*x + h_xyz*x*y) / 2.0;

  *det = 
    (*jacobian)[0]*((*jacobian)[4]*(*jacobian)[8] -
		    (*jacobian)[5]*(*jacobian)[7]) -
    (*jacobian)[1]*((*jacobian)[3]*(*jacobian)[8] -
		    (*jacobian)[5]*(*jacobian)[6]) +
    (*jacobian)[2]*((*jacobian)[3]*(*jacobian)[7] -
		    (*jacobian)[4]*(*jacobian)[6]);

  PetscLogFlops(152);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryHex3D::jacobian(double* jacobian,
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

  const double x4 = vertices[12];
  const double y4 = vertices[13];
  const double z4 = vertices[14];

  const double x5 = vertices[15];
  const double y5 = vertices[16];
  const double z5 = vertices[17];

  const double x6 = vertices[18];
  const double y6 = vertices[19];
  const double z6 = vertices[20];

  const double x7 = vertices[21];
  const double y7 = vertices[22];
  const double z7 = vertices[23];

  const double f_1 = (x1 - x0) / 2.0;
  const double g_1 = (y1 - y0) / 2.0;
  const double h_1 = (z1 - z0) / 2.0;

  const double f_3 = (x3 - x0) / 2.0;
  const double g_3 = (y3 - y0) / 2.0;
  const double h_3 = (z3 - z0) / 2.0;
  
  const double f_4 = (x4 - x0) / 2.0;
  const double g_4 = (y4 - y0) / 2.0;
  const double h_4 = (z4 - z0) / 2.0;

  const double f_01 = (x2 - x1 - x3 + x0) / 2.0;
  const double g_01 = (y2 - y1 - y3 + y0) / 2.0;
  const double h_01 = (z2 - z1 - z3 + z0) / 2.0;

  const double f_12 = (x7 - x3 - x4 + x0) / 2.0;
  const double g_12 = (y7 - y3 - y4 + y0) / 2.0;
  const double h_12 = (z7 - z3 - z4 + z0) / 2.0;

  const double f_02 = (x5 - x1 - x4 + x0) / 2.0;
  const double g_02 = (y5 - y1 - y4 + y0) / 2.0;
  const double h_02 = (z5 - z1 - z4 + z0) / 2.0;

  const double f_012 = (x6 - x0 + x1 - x2 + x3 + x4 - x5 - x7) / 2.0;
  const double g_012 = (y6 - y0 + y1 - y2 + y3 + y4 - y5 - y7) / 2.0;
  const double h_012 = (z6 - z0 + z1 - z2 + z3 + z4 - z5 - z7) / 2.0;

  for (int i=0, iR=0, iJ=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p1 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p2 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    assert(0 <= p2 && p2 <= 1.0);
    const double j0 = f_1 + f_01*p1 + f_02*p2 + f_012*p1*p2;
    const double j1 = f_3 + f_01*p0 + f_12*p2 + f_012*p0*p2;
    const double j2 = f_4 + f_12*p1 + f_02*p0 + f_012*p0*p1;
    
    const double j3 = g_1 + g_01*p1 + g_02*p2 + g_012*p1*p2;
    const double j4 = g_3 + g_01*p0 + g_12*p2 + g_012*p0*p2;
    const double j5 = g_4 + g_12*p1 + g_02*p0 + g_012*p0*p1;
    
    const double j6 = h_1 + h_01*p1 + h_02*p2 + h_012*p1*p2;
    const double j7 = h_3 + h_01*p0 + h_12*p2 + h_012*p0*p2;
    const double j8 = h_4 + h_12*p1 + h_02*p0 + h_012*p0*p1;

    jacobian[iJ++] = j0;
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    jacobian[iJ++] = j3;
    jacobian[iJ++] = j4;
    jacobian[iJ++] = j5;
    jacobian[iJ++] = j6;
    jacobian[iJ++] = j7;
    jacobian[iJ++] = j8;
    det[i] = j0*(j4*j8 - j5*j7) - j1*(j3*j8 - j5*j6) + j2*(j3*j7 - j4*j6);
  } // for

  PetscLogFlops(78 + npts*69);
} // jacobian


// End of file

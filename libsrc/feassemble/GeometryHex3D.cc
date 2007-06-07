// -*- C++ -*-
//
// ======================================================================
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ======================================================================
//

#include <portinfo>

#include "GeometryHex3D.hh" // implementation of class methods

#include "GeometryQuad3D.hh" // USES GeometryQuad3D

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryHex3D::GeometryHex3D(void) :
  CellGeometry(3, 3, 8)
{ // constructor
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
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryHex3D::jacobian(double_array* jacobian,
					  const double_array& vertices,
					  const double_array& location) const
{ // jacobian
  assert(0 != jacobian);

  assert(numCorners()*spaceDim() == vertices.size());
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

  const double x = location[0];
  const double y = location[1];
  const double z = location[2];
  assert(0 <= x && x <= 1.0);
  assert(0 <= y && y <= 1.0);
  assert(0 <= z && z <= 1.0);

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

  (*jacobian)[0] = x1 - x0 + f_xy*y + f_xz*z + f_xyz*y*z;
  (*jacobian)[1] = x3 - x0 + f_xy*x + f_yz*z + f_xyz*x*z;
  (*jacobian)[2] = x4 - x0 + f_yz*y + f_xz*x + f_xyz*x*y;

  (*jacobian)[3] = y1 - y0 + g_xy*y + g_xz*z + g_xyz*y*z;
  (*jacobian)[4] = y3 - y0 + g_xy*x + g_yz*z + g_xyz*x*z;
  (*jacobian)[5] = y4 - y0 + g_yz*y + g_xz*x + g_xyz*x*y;

  (*jacobian)[6] = z1 - z0 + h_xy*y + h_xz*z + h_xyz*y*z;
  (*jacobian)[7] = z3 - z0 + h_xy*x + h_yz*z + h_xyz*x*z;
  (*jacobian)[8] = z4 - z0 + h_yz*y + h_xz*x + h_xyz*x*y;
} // jacobian


// End of file

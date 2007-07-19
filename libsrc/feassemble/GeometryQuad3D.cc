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

#include "GeometryQuad3D.hh" // implementation of class methods

#include "GeometryLine3D.hh" // USES GeometryLine3D

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

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
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad3D::jacobian(double_array* jacobian,
					  double* det,
					  const double_array& vertices,
					  const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

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

  const double x = 0.5 * (location[0] + 1.0);
  const double y = 0.5 * (location[1] + 1.0);
  assert(0 <= x && x <= 1.0);
  assert(0 <= y && y <= 1.0);

  const double f_xy = x2 - x1 - x3 + x0;
  const double g_xy = y2 - y1 - y3 + y0;
  const double h_xy = z2 - z1 - z3 + z0;

  (*jacobian)[0] = x1 - x0 + f_xy*y;
  (*jacobian)[1] = x3 - x0 + f_xy*x;

  (*jacobian)[2] = y1 - y0 + g_xy*y;
  (*jacobian)[3] = y3 - y0 + g_xy*x;

  (*jacobian)[4] = z1 - z0 + h_xy*y;
  (*jacobian)[5] = z3 - z0 + h_xy*x;

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
} // jacobian


// End of file

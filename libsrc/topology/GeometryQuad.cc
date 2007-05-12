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

#include "GeometryQuad.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::GeometryQuad::GeometryQuad(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::topology::GeometryQuad::~GeometryQuad(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::topology::GeometryQuad::jacobian(double_array* jacobian,
					  const double_array& vertices,
					  const double_array& location) const
{ // jacobian
  assert(0 != jacobian);

  assert(2*4 == vertices.size());
  assert(2 == location.size());
  assert(4 == jacobian->size());
  
  const double x0 = vertices[0];
  const double y0 = vertices[1];

  const double x1 = vertices[2];
  const double y1 = vertices[3];

  const double x2 = vertices[4];
  const double y2 = vertices[5];

  const double x3 = vertices[6];
  const double y3 = vertices[7];

  const double x = location[0];
  const double y = location[1];

  const double f_xy = x2 - x1 - x3 + x0;
  const double g_xy = y2 - y1 - y3 + y0;

  (*jacobian)[0] = x1 - x0 + f_xy*y;
  (*jacobian)[1] = x3 - x0 + f_xy*x;
  (*jacobian)[2] = y1 - y0 + g_xy*y;
  (*jacobian)[3] = y3 - y0 + g_xy*x;
} // jacobian


// End of file

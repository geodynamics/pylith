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

  (*jacobian)[0] = x1 - x0 + (x3 - x1 - x2 + x0)*y;
  (*jacobian)[1] = x2 - x0 + (x3 - x1 - x2 + x0)*x;
  (*jacobian)[2] = y1 - y0 + (y3 - y1 - y2 + y0)*y;
  (*jacobian)[3] = y2 - y0 + (y3 - y1 - y2 + y0)*x;
} // jacobian


// End of file

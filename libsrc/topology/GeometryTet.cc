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

#include "GeometryTet.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::GeometryTet::GeometryTet(void)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::topology::GeometryTet::~GeometryTet(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::topology::GeometryTet::jacobian(double_array* jacobian,
					const double_array& vertices,
					const double_array& location) const
{ // jacobian
  assert(0 != jacobian);

  assert(12 == vertices.size());
  assert(3 == location.size());
  assert(9 == jacobian->size());
  
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

  (*jacobian)[0] = x1 - x0;
  (*jacobian)[1] = x2 - x0;
  (*jacobian)[2] = x3 - x0;
  (*jacobian)[3] = y1 - y0;
  (*jacobian)[4] = y2 - y0;
  (*jacobian)[5] = y3 - y0;
  (*jacobian)[6] = z1 - z0;
  (*jacobian)[7] = z2 - z0;
  (*jacobian)[8] = z3 - z0;
} // jacobian


// End of file

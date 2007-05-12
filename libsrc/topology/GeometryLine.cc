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

#include "GeometryLine.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::GeometryLine::GeometryLine(void) :
  CellGeometry(1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::topology::GeometryLine::~GeometryLine(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::topology::GeometryLine::jacobian(double_array* jacobian,
					  const double_array& vertices,
					  const double_array& location) const
{ // jacobian
  assert(0 != jacobian);

  assert(2 == vertices.size());
  assert(1 == location.size());
  assert(1 == jacobian->size());

  (*jacobian)[0] = vertices[1]-vertices[0];
} // jacobian


// End of file

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

#include "GeometryLine2D.hh" // implementation of class methods

#include "GeometryPoint2D.hh" // USES GeometryPoint

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::GeometryLine2D::GeometryLine2D(void) :
  CellGeometry(1, 2, 2)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::topology::GeometryLine2D::~GeometryLine2D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::topology::CellGeometry*
pylith::topology::GeometryLine2D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryPoint2D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::topology::GeometryLine2D::jacobian(double_array* jacobian,
					   const double_array& vertices,
					   const double_array& location) const
{ // jacobian
  assert(0 != jacobian);

  assert(numCorners()*spaceDim() == vertices.size());
  assert(spaceDim()*cellDim() == jacobian->size());

  const double x0 = vertices[0];
  const double y0 = vertices[1];

  const double x1 = vertices[2];
  const double y1 = vertices[3];

  (*jacobian)[0] = x1 - x0;
  (*jacobian)[1] = y1 - y0;
} // jacobian


// End of file

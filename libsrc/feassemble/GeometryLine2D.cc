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
pylith::feassemble::GeometryLine2D::GeometryLine2D(void) :
  CellGeometry(LINE, 2)
{ // constructor
  const double vertices[] = {
    -1.0,
    +1.0,
  };
  _setVertices(vertices, 2, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryLine2D::~GeometryLine2D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine2D::clone(void) const
{ // clone
  return new GeometryLine2D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine2D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryPoint2D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine2D::jacobian(double_array* jacobian,
					   double* det,
					   const double_array& vertices,
					   const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert(numCorners()*spaceDim() == vertices.size());
  assert(spaceDim()*cellDim() == jacobian->size());

  const double x0 = vertices[0];
  const double y0 = vertices[1];

  const double x1 = vertices[2];
  const double y1 = vertices[3];

  (*jacobian)[0] = x1 - x0;
  (*jacobian)[1] = y1 - y0;
  *det = sqrt(pow((*jacobian)[0], 2) +
	      pow((*jacobian)[1], 2));
} // jacobian


// End of file

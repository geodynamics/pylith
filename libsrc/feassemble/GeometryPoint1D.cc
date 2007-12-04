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

#include "GeometryPoint1D.hh" // implementation of class methods

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryPoint1D::GeometryPoint1D(void) :
  CellGeometry(POINT, 1)
{ // constructor
  const double vertices[] = { 0.0 };
  _setVertices(vertices, 1, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryPoint1D::~GeometryPoint1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryPoint1D::clone(void) const
{ // clone
  return new GeometryPoint1D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryPoint1D::geometryLowerDim(void) const
{ // geometryLowerDim
  return 0;
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryPoint1D::coordsRefToGlobal(double* coordsGlobal,
						       const double* coordsRef,
						       const double* vertices,
						       const int dim) const
{ // coordsRefToGlobal
  assert(0 != coordsGlobal);
  assert(0 != coordsRef);
  assert(0 != vertices);
  assert(1 == dim);
  assert(spaceDim() == dim);

  coordsGlobal[0] = vertices[0];
} // coordsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryPoint1D::jacobian(double_array* jacobian,
					    double* det,
					    const double_array& vertices,
					    const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert(1 == jacobian->size());
  
  (*jacobian)[0] = 1.0;
  *det = 1.0;
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryPoint1D::jacobian(double* jacobian,
					      double* det,
					      const double* vertices,
					      const double* location,
					      const int dim) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != location);
  assert(1 == dim);
  assert(spaceDim() == dim);

  jacobian[0] = 1.0;
  *det = 1.0;
} // jacobian


// End of file

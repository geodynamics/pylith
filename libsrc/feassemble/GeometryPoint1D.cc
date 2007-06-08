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
  CellGeometry(0, 1, 1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryPoint1D::~GeometryPoint1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::GeometryPoint1D::GeometryPoint1D(const GeometryPoint1D& g) :
  CellGeometry(0, 1, 1)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryPoint1D::clone(void) const
{ // clone
  return new GeometryPoint1D(*this);
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryPoint1D::geometryLowerDim(void) const
{ // geometryLowerDim
  return 0;
} // geometryLowerDim

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


// End of file

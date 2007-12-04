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

#include "GeometryLine1D.hh" // implementation of class methods

#include "GeometryPoint1D.hh" // USES GeometryPoint

#include "pylith/utils/array.hh" // USES double_array

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryLine1D::GeometryLine1D(void) :
  CellGeometry(LINE, 1)
{ // constructor
  const double vertices[] = {
    -1.0,
    +1.0,
  };
  _setVertices(vertices, 2, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryLine1D::~GeometryLine1D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine1D::clone(void) const
{ // clone
  return new GeometryLine1D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine1D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryPoint1D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryLine1D::coordsRefToGlobal(double* coordsGlobal,
						      const double* coordsRef,
						      const double* vertices,
						      const int dim) const
{ // coordsRefToGlobal
  assert(0 != coordsGlobal);
  assert(0 != coordsRef);
  assert(0 != vertices);
  assert(1 == dim);
  assert(spaceDim() == dim);

  const double x0 = vertices[0];
  const double x1 = vertices[1];

  const double p0 = 0.5*(1.0+coordsRef[0]);
  coordsGlobal[0] = x0 + (x1-x0) * p0;

  PetscLogFlopsNoCheck(5);
} // coordsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine1D::jacobian(double_array* jacobian,
					     double* det,
					     const double_array& vertices,
					     const double_array& location) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);

  assert(numCorners()*spaceDim() == vertices.size());
  assert(spaceDim()*cellDim() == jacobian->size());

  const double x0 = vertices[0];
  const double x1 = vertices[1];

  (*jacobian)[0] = (x1 - x0)/2.0;
  *det = (*jacobian)[0];
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine1D::jacobian(double* jacobian,
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

  const double x0 = vertices[0];
  const double x1 = vertices[1];

  jacobian[0] = (x1 - x0)/2.0;
  *det = jacobian[0];
  PetscLogFlopsNoCheck(2);
} // jacobian


// End of file

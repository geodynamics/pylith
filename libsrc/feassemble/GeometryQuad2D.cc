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

#include "GeometryQuad2D.hh" // implementation of class methods

#include "GeometryLine2D.hh" // USES GeometryLine2D

#include "petsc.h" // USES PetscLogFlopsNoCheck

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryQuad2D::GeometryQuad2D(void) :
  CellGeometry(QUADRILATERAL, 2)
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
pylith::feassemble::GeometryQuad2D::~GeometryQuad2D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryQuad2D::clone(void) const
{ // clone
  return new GeometryQuad2D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryQuad2D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryLine2D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryQuad2D::coordsRefToGlobal(double* coordsGlobal,
						      const double* coordsRef,
						      const double* vertices,
						      const int dim) const
{ // coordsRefToGlobal
  assert(0 != coordsGlobal);
  assert(0 != coordsRef);
  assert(0 != vertices);
  assert(2 == dim);
  assert(spaceDim() == dim);

  const double x0 = vertices[0];
  const double y0 = vertices[1];

  const double x1 = vertices[2];
  const double y1 = vertices[3];

  const double x2 = vertices[4];
  const double y2 = vertices[5];

  const double x3 = vertices[6];
  const double y3 = vertices[7];

  const double p0 = 0.5*(1.0+coordsRef[0]);
  const double p1 = 0.5*(1.0+coordsRef[1]);
  const double f_01 = x2 - x1 - x3 + x0;
  const double g_01 = y2 - y1 - y3 + y0;
  coordsGlobal[0] = x0 + (x1-x0) * p0 + (x2-x0) * p1 + f_01 * p0 * p1;
  coordsGlobal[1] = y0 + (y1-y0) * p0 + (y2-y0) * p1 + g_01 * p0 * p1;

  PetscLogFlopsNoCheck(28);
} // coordsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad2D::jacobian(double_array* jacobian,
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

  const double x1 = vertices[2];
  const double y1 = vertices[3];

  const double x2 = vertices[4];
  const double y2 = vertices[5];

  const double x3 = vertices[6];
  const double y3 = vertices[7];

  const double x = 0.5 * (location[0] + 1.0);
  const double y = 0.5 * (location[1] + 1.0);
  assert(0 <= x && x <= 1.0);
  assert(0 <= y && y <= 1.0);

  const double f_xy = x2 - x1 - x3 + x0;
  const double g_xy = y2 - y1 - y3 + y0;

  (*jacobian)[0] = (x1 - x0 + f_xy*y) / 2.0;
  (*jacobian)[1] = (x3 - x0 + f_xy*x) / 2.0;
  (*jacobian)[2] = (y1 - y0 + g_xy*y) / 2.0;
  (*jacobian)[3] = (y3 - y0 + g_xy*x) / 2.0;

  *det = 
    (*jacobian)[0]*(*jacobian)[3] - 
    (*jacobian)[1]*(*jacobian)[2];

  PetscLogFlopsNoCheck(31);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad2D::jacobian(double* jacobian,
					     double* det,
					     const double* vertices,
					     const double* location,
					     const int dim) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != location);
  assert(2 == dim);
  assert(spaceDim() == dim);
    
  const double x0 = vertices[0];
  const double y0 = vertices[1];

  const double x1 = vertices[2];
  const double y1 = vertices[3];

  const double x2 = vertices[4];
  const double y2 = vertices[5];

  const double x3 = vertices[6];
  const double y3 = vertices[7];

  const double x = 0.5 * (location[0] + 1.0);
  const double y = 0.5 * (location[1] + 1.0);
  assert(0 <= x && x <= 1.0);
  assert(0 <= y && y <= 1.0);

  const double f_xy = x2 - x1 - x3 + x0;
  const double g_xy = y2 - y1 - y3 + y0;

  jacobian[0] = (x1 - x0 + f_xy*y) / 2.0;
  jacobian[1] = (x3 - x0 + f_xy*x) / 2.0;
  jacobian[2] = (y1 - y0 + g_xy*y) / 2.0;
  jacobian[3] = (y3 - y0 + g_xy*x) / 2.0;

  *det = 
    jacobian[0]*jacobian[3] - 
    jacobian[1]*jacobian[2];
  PetscLogFlopsNoCheck(29);
} // jacobian


// End of file

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

#include "GeometryTet3D.hh" // implementation of class methods

#include "GeometryTri3D.hh" // USES GeometryTri3D

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryTet3D::GeometryTet3D(void) :
  CellGeometry(3, 3, 4)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryTet3D::~GeometryTet3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::GeometryTet3D::GeometryTet3D(const GeometryTet3D& g) :
  CellGeometry(3, 3, 4)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTet3D::clone(void) const
{ // clone
  return new GeometryTet3D(*this);
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTet3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryTri3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTet3D::jacobian(double_array* jacobian,
					  const double_array& vertices,
					  const double_array& location) const
{ // jacobian
  assert(0 != jacobian);

  assert(numCorners()*spaceDim() == vertices.size());
  assert(spaceDim()*cellDim() == jacobian->size());
  
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

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

#include "GeometryTri3D.hh" // implementation of class methods

#include "GeometryLine3D.hh" // USES GeometryLine3D

#include "pylith/utils/array.hh" // USES double_array

#include <assert.h> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::topology::GeometryTri3D::GeometryTri3D(void) :
  CellGeometry(2, 3, 3)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::topology::GeometryTri3D::~GeometryTri3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::topology::CellGeometry*
pylith::topology::GeometryTri3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryLine3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::topology::GeometryTri3D::jacobian(double_array* jacobian,
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

  (*jacobian)[0] = x1 - x0;
  (*jacobian)[1] = x2 - x0;

  (*jacobian)[2] = y1 - y0;
  (*jacobian)[3] = y2 - y0;

  (*jacobian)[4] = z1 - z0;
  (*jacobian)[5] = z2 - z0;
} // jacobian


// End of file

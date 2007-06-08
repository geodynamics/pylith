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
pylith::feassemble::GeometryTri3D::GeometryTri3D(void) :
  CellGeometry(2, 3, 3)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryTri3D::~GeometryTri3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::GeometryTri3D::GeometryTri3D(const GeometryTri3D& g) :
  CellGeometry(2, 3, 3)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTri3D::clone(void) const
{ // clone
  return new GeometryTri3D(*this);
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTri3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryLine3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTri3D::jacobian(double_array* jacobian,
					  double* det,
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

  const double jj00 = 
    (*jacobian)[0]*(*jacobian)[0] +
    (*jacobian)[2]*(*jacobian)[2] +
    (*jacobian)[4]*(*jacobian)[4];
  const double jj10 =
    (*jacobian)[0]*(*jacobian)[1] +
    (*jacobian)[2]*(*jacobian)[3] +
    (*jacobian)[4]*(*jacobian)[5];
  const double jj01 = jj10;
  const double jj11 = 
    (*jacobian)[1]*(*jacobian)[1] +
    (*jacobian)[3]*(*jacobian)[3] +
    (*jacobian)[5]*(*jacobian)[5];
  *det = sqrt(jj00*jj11 - jj01*jj10);
} // jacobian


// End of file

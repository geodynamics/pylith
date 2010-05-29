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

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES double_array

#include <cassert> // USES assert()

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
pylith::feassemble::GeometryQuad2D::ptsRefToGlobal(double* ptsGlobal,
						   const double* ptsRef,
						   const double* vertices,
						   const int dim,
						   const int npts) const
{ // ptsRefToGlobal
  assert(0 != ptsGlobal);
  assert(0 != ptsRef);
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

  const double f_1 = x1 - x0;
  const double g_1 = y1 - y0;

  const double f_3 = x3 - x0;
  const double g_3 = y3 - y0;

  const double f_01 = x2 - x1 - x3 + x0;
  const double g_01 = y2 - y1 - y3 + y0;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p1 = 0.5 * (1.0 + ptsRef[iR++]);
    ptsGlobal[iG++] = x0 + f_1 * p0 + f_3 * p1 + f_01 * p0 * p1;
    ptsGlobal[iG++] = y0 + g_1 * p0 + g_3 * p1 + g_01 * p0 * p1;
  } // for

  PetscLogFlops(10 + npts*18);
} // ptsRefToGlobal

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

  PetscLogFlops(31);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryQuad2D::jacobian(double* jacobian,
					     double* det,
					     const double* vertices,
					     const double* ptsRef,
					     const int dim,
					     const int npts) const
{ // jacobian
  assert(0 != jacobian);
  assert(0 != det);
  assert(0 != vertices);
  assert(0 != ptsRef);
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

  const double f_1 = x1 - x0;
  const double g_1 = y1 - y0;

  const double f_3 = x3 - x0;
  const double g_3 = y3 - y0;

  const double f_01 = x2 - x1 - x3 + x0;
  const double g_01 = y2 - y1 - y3 + y0;

  for (int i=0, iR=0, iJ=0; i < npts; ++i) {
    const double p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const double p1 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    const double j00 = (f_1 + f_01 * p1) / 2.0; 
    const double j01 = (f_3 + f_01 * p0) / 2.0; 
    const double j10 = (g_1 + g_01 * p1) / 2.0;
    const double j11 = (g_3 + g_01 * p0) / 2.0; 

    jacobian[iJ++] = j00;
    jacobian[iJ++] = j01;
    jacobian[iJ++] = j10;
    jacobian[iJ++] = j11;
    det[i] = j00*j11 - j01*j10;
  } // for

  PetscLogFlops(10 + npts*19);
} // jacobian


// End of file

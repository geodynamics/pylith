// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "GeometryTet3D.hh" // implementation of class methods

#include "GeometryTri3D.hh" // USES GeometryTri3D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES scalar_array
#include "pylith/utils/error.h" // USES std::logic_error

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryTet3D::GeometryTet3D(void) :
  CellGeometry(TETRAHEDRON, 3)
{ // constructor
  const PylithScalar vertices[] = {
    -1.0,  -1.0,  -1.0,
    +1.0,  -1.0,  -1.0,
    -1.0,  +1.0,  -1.0,
    -1.0,  -1.0,  +1.0,
  };
  _setVertices(vertices, 4, 3);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryTet3D::~GeometryTet3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTet3D::clone(void) const
{ // clone
  return new GeometryTet3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTet3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryTri3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryTet3D::ptsRefToGlobal(PylithScalar* ptsGlobal,
						  const PylithScalar* ptsRef,
						  const PylithScalar* vertices,
						  const int dim,
						  const int npts) const
{ // ptsRefToGlobal
  assert(ptsGlobal);
  assert(ptsRef);
  assert(vertices);
  assert(3 == dim);
  assert(spaceDim() == dim);

  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];
  const PylithScalar z0 = vertices[2];

  const PylithScalar x1 = vertices[3];
  const PylithScalar y1 = vertices[4];
  const PylithScalar z1 = vertices[5];

  const PylithScalar x2 = vertices[6];
  const PylithScalar y2 = vertices[7];
  const PylithScalar z2 = vertices[8];

  const PylithScalar x3 = vertices[9];
  const PylithScalar y3 = vertices[10];
  const PylithScalar z3 = vertices[11];

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;
  const PylithScalar h_1 = z1 - z0;

  const PylithScalar f_2 = x2 - x0;
  const PylithScalar g_2 = y2 - y0;
  const PylithScalar h_2 = z2 - z0;

  const PylithScalar f_3 = x3 - x0;
  const PylithScalar g_3 = y3 - y0;
  const PylithScalar h_3 = z3 - z0;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p1 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p2 = 0.5 * (1.0 + ptsRef[iR++]);
    assert(0 <= p0 && p0 <= 1.0);
    assert(0 <= p1 && p1 <= 1.0);
    assert(0 <= p2 && p2 <= 1.0);

    ptsGlobal[iG++] = x0 + f_1 * p0 + f_2 * p1 + f_3 * p2;
    ptsGlobal[iG++] = y0 + g_1 * p0 + g_2 * p1 + g_3 * p2;;
    ptsGlobal[iG++] = z0 + h_1 * p0 + h_2 * p1 + h_3 * p2;
  } // for

  PetscLogFlops(9 + npts*24);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTet3D::jacobian(scalar_array* jacobian,
					    PylithScalar* det,
					    const PylithScalar* vertices,
					    const int numVertices,
					    const int spaceDim,
					    const PylithScalar* location,
					    const int cellDim) const
{ // jacobian
  assert(jacobian);
  assert(det);
  assert(vertices);
  assert(location);

  assert(this->numCorners() == numVertices || // linear
	 this->numCorners()+1 == numVertices); // quadratic
  assert(this->spaceDim() == spaceDim);
  assert(this->cellDim() == cellDim);
  assert(size_t(spaceDim*cellDim) == jacobian->size());
  
  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];
  const PylithScalar z0 = vertices[2];

  const PylithScalar x2 = vertices[3];
  const PylithScalar y2 = vertices[4];
  const PylithScalar z2 = vertices[5];

  const PylithScalar x1 = vertices[6];
  const PylithScalar y1 = vertices[7];
  const PylithScalar z1 = vertices[8];

  const PylithScalar x3 = vertices[9];
  const PylithScalar y3 = vertices[10];
  const PylithScalar z3 = vertices[11];

  (*jacobian)[0] = (x1 - x0) / 2.0;
  (*jacobian)[1] = (x2 - x0) / 2.0;
  (*jacobian)[2] = (x3 - x0) / 2.0;
  (*jacobian)[3] = (y1 - y0) / 2.0;
  (*jacobian)[4] = (y2 - y0) / 2.0;
  (*jacobian)[5] = (y3 - y0) / 2.0;
  (*jacobian)[6] = (z1 - z0) / 2.0;
  (*jacobian)[7] = (z2 - z0) / 2.0;
  (*jacobian)[8] = (z3 - z0) / 2.0;

  *det = 
    (*jacobian)[0]*((*jacobian)[4]*(*jacobian)[8] -
		    (*jacobian)[5]*(*jacobian)[7]) -
    (*jacobian)[1]*((*jacobian)[3]*(*jacobian)[8] -
		    (*jacobian)[5]*(*jacobian)[6]) +
    (*jacobian)[2]*((*jacobian)[3]*(*jacobian)[7] -
		    (*jacobian)[4]*(*jacobian)[6]);

  PetscLogFlops(32);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTet3D::jacobian(PylithScalar* jacobian,
					    PylithScalar* det,
					    const PylithScalar* vertices,
					    const PylithScalar* ptsRef,
					    const int dim,
					    const int npts) const
{ // jacobian
  assert(jacobian);
  assert(det);
  assert(vertices);
  assert(ptsRef);
  assert(3 == dim);
  assert(spaceDim() == dim);
  
  const PylithScalar x0 = vertices[0];
  const PylithScalar y0 = vertices[1];
  const PylithScalar z0 = vertices[2];

  const PylithScalar x1 = vertices[3];
  const PylithScalar y1 = vertices[4];
  const PylithScalar z1 = vertices[5];

  const PylithScalar x2 = vertices[6];
  const PylithScalar y2 = vertices[7];
  const PylithScalar z2 = vertices[8];

  const PylithScalar x3 = vertices[9];
  const PylithScalar y3 = vertices[10];
  const PylithScalar z3 = vertices[11];

  const PylithScalar j0 = (x1 - x0) / 2.0;
  const PylithScalar j1 = (x2 - x0) / 2.0;
  const PylithScalar j2 = (x3 - x0) / 2.0;

  const PylithScalar j3 = (y1 - y0) / 2.0;
  const PylithScalar j4 = (y2 - y0) / 2.0;
  const PylithScalar j5 = (y3 - y0) / 2.0;

  const PylithScalar j6 = (z1 - z0) / 2.0;
  const PylithScalar j7 = (z2 - z0) / 2.0;
  const PylithScalar j8 = (z3 - z0) / 2.0;

  const PylithScalar jdet = 
    j0*(j4*j8 - j5*j7) -
    j1*(j3*j8 - j5*j6) +
    j2*(j3*j7 - j4*j6);

  for (int i=0, iJ=0; i < npts; ++i) {
    jacobian[iJ++] = j0;
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    jacobian[iJ++] = j3;
    jacobian[iJ++] = j4;
    jacobian[iJ++] = j5;
    jacobian[iJ++] = j6;
    jacobian[iJ++] = j7;
    jacobian[iJ++] = j8;
    det[i] = jdet;
  } // for

  PetscLogFlops(32);
} // jacobian


// ----------------------------------------------------------------------
// Compute minimum width across cell.
PylithScalar
pylith::feassemble::GeometryTet3D::minCellWidth(const PylithScalar* coordinatesCell,
						const int numVertices,
						const int spaceDim) const
{ // minCellWidth
  const int numCorners = 4;
  const int dim = 3;

  assert(numCorners == numVertices || // linear
	 numCorners+6 == numVertices); // quadratic
  assert(dim == spaceDim);

  const int numEdges = 6;
  const int edges[numEdges][2] = {
    {0, 1}, {1, 2}, {2, 0},
    {0, 3}, {1, 3}, {2, 3},
  };

  PylithScalar minWidth = PYLITH_MAXSCALAR;
  for (int iedge=0; iedge < numEdges; ++iedge) {
    const int iA = edges[iedge][0];
    const int iB = edges[iedge][1];
    const PylithScalar xA = coordinatesCell[dim*iA  ];
    const PylithScalar yA = coordinatesCell[dim*iA+1];
    const PylithScalar zA = coordinatesCell[dim*iA+2];
    const PylithScalar xB = coordinatesCell[dim*iB  ];
    const PylithScalar yB = coordinatesCell[dim*iB+1];
    const PylithScalar zB = coordinatesCell[dim*iB+2];
    
    const PylithScalar edgeLen = sqrt(pow(xB-xA,2) + pow(yB-yA,2) + pow(zB-zA,2));
    if (edgeLen < minWidth) {
      minWidth = edgeLen;
    } // if
  } // for

  PetscLogFlops(numEdges*9);

  // Radius of inscribed sphere
  const PylithScalar v = volume(coordinatesCell, numVertices, spaceDim);
  const PylithScalar a = faceArea(coordinatesCell, numVertices, spaceDim, 0) +
    faceArea(coordinatesCell, numVertices, spaceDim, 1) +
    faceArea(coordinatesCell, numVertices, spaceDim, 2) +
    faceArea(coordinatesCell, numVertices, spaceDim, 3);
    
  const PylithScalar r = 3.0 * v / a;
  const PylithScalar rwidth = 6.38*r; // based on empirical tests
  if (rwidth < minWidth) {
    minWidth = rwidth;
  } // if

  PetscLogFlops(3);

  return minWidth;
} // minCellWidth

// ----------------------------------------------------------------------
// Compute cell volume.
PylithScalar
pylith::feassemble::GeometryTet3D::volume(const PylithScalar* coordinatesCell,
					  const int numVertices,
					  const int spaceDim) const
{ // volume
  const int numCorners = 4;
  const int dim = 3;

  assert(numCorners == numVertices || // linear
	 numCorners+6 == numVertices); // quadratic
  assert(dim == spaceDim);
  
  const PylithScalar x0 = coordinatesCell[0];
  const PylithScalar y0 = coordinatesCell[1];
  const PylithScalar z0 = coordinatesCell[2];
  
  const PylithScalar x2 = coordinatesCell[3];
  const PylithScalar y2 = coordinatesCell[4];
  const PylithScalar z2 = coordinatesCell[5];
  
  const PylithScalar x1 = coordinatesCell[6];
  const PylithScalar y1 = coordinatesCell[7];
  const PylithScalar z1 = coordinatesCell[8];
  
  const PylithScalar x3 = coordinatesCell[9];
  const PylithScalar y3 = coordinatesCell[10];
  const PylithScalar z3 = coordinatesCell[11];

  const PylithScalar det = 
    x1*(y2*z3-y3*z2)-y1*(x2*z3-x3*z2)+(x2*y3-x3*y2)*z1 - 
    x0*((y2*z3-y3*z2)-y1*(z3-z2)+(y3-y2)*z1) +
    y0*((x2*z3-x3*z2)-x1*(z3-z2)+(x3-x2)*z1) -
    z0*((x2*y3-x3*y2)-x1*(y3-y2)+(x3-x2)*y1);
  assert(det > 0.0);

  const PylithScalar v = det / 6.0;
  PetscLogFlops(48);
  
  return v;  
} // volume

// ----------------------------------------------------------------------
// Compute area of face.
PylithScalar
pylith::feassemble::GeometryTet3D::faceArea(const PylithScalar* coordinatesCell,
					    const int numVertices,
					    const int spaceDim,
					    const int face) const
{ // faceArea
  const int numCorners = 4;
  const int dim = 3;

  assert(numCorners == numVertices || // linear
	 numCorners+6 == numVertices); // quadratic
  assert(dim == spaceDim);

  const PylithScalar x0 = coordinatesCell[0];
  const PylithScalar y0 = coordinatesCell[1];
  const PylithScalar z0 = coordinatesCell[2];
  
  const PylithScalar x1 = coordinatesCell[3];
  const PylithScalar y1 = coordinatesCell[4];
  const PylithScalar z1 = coordinatesCell[5];
  
  const PylithScalar x2 = coordinatesCell[6];
  const PylithScalar y2 = coordinatesCell[7];
  const PylithScalar z2 = coordinatesCell[8];
  
  const PylithScalar x3 = coordinatesCell[9];
  const PylithScalar y3 = coordinatesCell[10];
  const PylithScalar z3 = coordinatesCell[11];

  PylithScalar a[3];
  PylithScalar b[3];
  switch (face) {
  case 0:
    a[0] = x3-x1;
    a[1] = y3-y1;
    a[2] = z3-z1;
    b[0] = x2-x1;
    b[1] = y2-y1;
    b[2] = z2-z1;
    break;
  case 1:
    a[0] = x3-x2;
    a[1] = y3-y2;
    a[2] = z3-z2;
    b[0] = x0-x2;
    b[1] = y0-y2;
    b[2] = z0-z2;
    break;
  case 2:
    a[0] = x3-x0;
    a[1] = y3-y0;
    a[2] = z3-z0;
    b[0] = x1-x0;
    b[1] = y1-y0;
    b[2] = z1-z0;
    break;
  case 3:
    a[0] = x1-x0;
    a[1] = y1-y0;
    a[2] = z1-z0;
    b[0] = x2-x0;
    b[1] = y2-y0;
    b[2] = z2-z0;
    break;
  default:
    assert(0);
    throw std::logic_error("Unknown face.");
  } // switch

  const PylithScalar areaX = a[1]*b[2] - a[2]*b[1];
  const PylithScalar areaY = a[2]*b[0] - a[0]*b[2];
  const PylithScalar areaZ = a[0]*b[1] - a[1]*b[0];

  const PylithScalar area = 0.5*sqrt(areaX*areaX + areaY*areaY + areaZ*areaZ);
  PetscLogFlops(22);
  
  return area;
} // faceArea


// End of file

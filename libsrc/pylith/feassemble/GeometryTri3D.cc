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

#include "GeometryTri3D.hh" // implementation of class methods

#include "GeometryLine3D.hh" // USES GeometryLine3D

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES scalar_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryTri3D::GeometryTri3D(void) :
  CellGeometry(TRIANGLE, 3)
{ // constructor
  const PylithScalar vertices[] = {
    -1.0,  -1.0,
    +1.0,  -1.0,
    -1.0,  +1.0,
  };
  _setVertices(vertices, 3, 2);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryTri3D::~GeometryTri3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTri3D::clone(void) const
{ // clone
  return new GeometryTri3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryTri3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return new GeometryLine3D();
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryTri3D::ptsRefToGlobal(PylithScalar* ptsGlobal,
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

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;
  const PylithScalar h_1 = z1 - z0;

  const PylithScalar f_2 = x2 - x0;
  const PylithScalar g_2 = y2 - y0;
  const PylithScalar h_2 = z2 - z0;

  for (int i=0, iR=0, iG=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[iR++]);
    const PylithScalar p1 = 0.5 * (1.0 + ptsRef[iR++]);
    ptsGlobal[iG++] = x0 + f_1 * p0 + f_2 * p1;
    ptsGlobal[iG++] = y0 + g_1 * p0 + g_2 * p1;
    ptsGlobal[iG++] = z0 + h_1 * p0 + h_2 * p1;
  } // for

  PetscLogFlops(22);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTri3D::jacobian(scalar_array* jacobian,
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

  const PylithScalar x1 = vertices[3];
  const PylithScalar y1 = vertices[4];
  const PylithScalar z1 = vertices[5];

  const PylithScalar x2 = vertices[6];
  const PylithScalar y2 = vertices[7];
  const PylithScalar z2 = vertices[8];

  (*jacobian)[0] = (x1 - x0) / 2.0;
  (*jacobian)[1] = (x2 - x0) / 2.0;

  (*jacobian)[2] = (y1 - y0) / 2.0;
  (*jacobian)[3] = (y2 - y0) / 2.0;

  (*jacobian)[4] = (z1 - z0) / 2.0;
  (*jacobian)[5] = (z2 - z0) / 2.0;

  const PylithScalar jj00 = 
    (*jacobian)[0]*(*jacobian)[0] +
    (*jacobian)[2]*(*jacobian)[2] +
    (*jacobian)[4]*(*jacobian)[4];
  const PylithScalar jj10 =
    (*jacobian)[0]*(*jacobian)[1] +
    (*jacobian)[2]*(*jacobian)[3] +
    (*jacobian)[4]*(*jacobian)[5];
  const PylithScalar jj01 = jj10;
  const PylithScalar jj11 = 
    (*jacobian)[1]*(*jacobian)[1] +
    (*jacobian)[3]*(*jacobian)[3] +
    (*jacobian)[5]*(*jacobian)[5];
  *det = sqrt(jj00*jj11 - jj01*jj10);
  PetscLogFlops(25);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryTri3D::jacobian(PylithScalar* jacobian,
					    PylithScalar* det,
					    const PylithScalar* vertices,
					    const PylithScalar* location,
					    const int dim,
					    const int npts) const
{ // jacobian
  assert(jacobian);
  assert(det);
  assert(vertices);
  assert(location);
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

  const PylithScalar j00 = (x1 - x0) / 2.0;
  const PylithScalar j01 = (x2 - x0) / 2.0;

  const PylithScalar j10 = (y1 - y0) / 2.0;
  const PylithScalar j11 = (y2 - y0) / 2.0;

  const PylithScalar j20 = (z1 - z0) / 2.0;
  const PylithScalar j21 = (z2 - z0) / 2.0;

  const PylithScalar jj00 = j00*j00 + j10*j10 + j20*j20;
  const PylithScalar jj10 = j00*j01 + j10*j11 + j20*j21;
  const PylithScalar jj01 = jj10;
  const PylithScalar jj11 = j01*j01 + j11*j11 + j21*j21;
  const PylithScalar jdet = sqrt(jj00*jj11 - jj01*jj10);

  for (int i=0, iJ=0; i < npts; ++i) {
    jacobian[iJ++] = j00;
    jacobian[iJ++] = j01;
    jacobian[iJ++] = j10;
    jacobian[iJ++] = j11;
    jacobian[iJ++] = j20;
    jacobian[iJ++] = j21;
    det[i] = jdet;
  } // for

  PetscLogFlops(31);
} // jacobian


// ----------------------------------------------------------------------
// Compute minimum width across cell.
PylithScalar
pylith::feassemble::GeometryTri3D::minCellWidth(const PylithScalar* coordinatesCell,
						const int numVertices,
						const int spaceDim) const
{ // minCellWidth
  const int numCorners = 2;
  const int dim = 3;

  assert(numCorners == numVertices);
  assert(dim == spaceDim);

  const int numEdges = 3;
  const int edges[numEdges][2] = {
    {0, 1}, {1, 2}, {2, 0},
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

#if 1
  // Ad-hoc to account for distorted cells.
  // Radius of inscribed circle.
  const PylithScalar xA = coordinatesCell[0];
  const PylithScalar yA = coordinatesCell[1];
  const PylithScalar zA = coordinatesCell[2];
  const PylithScalar xB = coordinatesCell[3];
  const PylithScalar yB = coordinatesCell[4];
  const PylithScalar zB = coordinatesCell[5];
  const PylithScalar xC = coordinatesCell[6];
  const PylithScalar yC = coordinatesCell[7];
  const PylithScalar zC = coordinatesCell[8];

  const PylithScalar c  = sqrt(pow(xB-xA,2) + pow(yB-yA, 2) + pow(zB-zA, 2));
  const PylithScalar a  = sqrt(pow(xC-xB,2) + pow(yC-yB, 2) + pow(zC-zB, 2));
  const PylithScalar b  = sqrt(pow(xA-xC,2) + pow(yA-yC, 2) + pow(zA-zC, 2));
  
  const PylithScalar k = (a + b + c) / 3.0;
  const PylithScalar r = sqrt(k*(k-a)*(k-b)*(k-c)) / k;
  minWidth = r;

#endif

  return minWidth;
} // minCellWidth


// End of file

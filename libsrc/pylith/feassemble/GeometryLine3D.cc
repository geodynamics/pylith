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

#include "GeometryLine3D.hh" // implementation of class methods

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/constdefs.h" // USES scalar_array

#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::GeometryLine3D::GeometryLine3D(void) :
  CellGeometry(LINE, 3)
{ // constructor
  const PylithScalar vertices[] = {
    -1.0,
    +1.0,
  };
  _setVertices(vertices, 2, 1);
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::GeometryLine3D::~GeometryLine3D(void)
{ // destructor
} // destructor

// ----------------------------------------------------------------------
// Create a copy of geometry.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine3D::clone(void) const
{ // clone
  return new GeometryLine3D();
} // clone

// ----------------------------------------------------------------------
// Get cell geometry for lower dimension cell.
pylith::feassemble::CellGeometry*
pylith::feassemble::GeometryLine3D::geometryLowerDim(void) const
{ // geometryLowerDim
  return 0;
} // geometryLowerDim

// ----------------------------------------------------------------------
// Transform coordinates in reference cell to global coordinates.
void
pylith::feassemble::GeometryLine3D::ptsRefToGlobal(PylithScalar* ptsGlobal,
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

  const PylithScalar f_1 = x1 - x0;
  const PylithScalar g_1 = y1 - y0;
  const PylithScalar h_1 = z1 - z0;

  for (int i=0, iG=0; i < npts; ++i) {
    const PylithScalar p0 = 0.5 * (1.0 + ptsRef[i]);
    ptsGlobal[iG++] = x0 + f_1 * p0;
    ptsGlobal[iG++] = y0 + g_1 * p0;
    ptsGlobal[iG++] = z0 + h_1 * p0;
  } // for

  PetscLogFlops(3 + npts*8);
} // ptsRefToGlobal

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine3D::jacobian(scalar_array* jacobian,
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

  (*jacobian)[0] = (x1 - x0)/2.0;
  (*jacobian)[1] = (y1 - y0)/2.0;
  (*jacobian)[2] = (z1 - z0)/2.0;
  *det = sqrt(pow((*jacobian)[0], 2) +
	      pow((*jacobian)[1], 2) +
	      pow((*jacobian)[2], 2));
  PetscLogFlops(12);
} // jacobian

// ----------------------------------------------------------------------
// Compute Jacobian at location in cell.
void
pylith::feassemble::GeometryLine3D::jacobian(PylithScalar* jacobian,
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

  const PylithScalar j1 = (x1 - x0) / 2.0;
  const PylithScalar j2 = (y1 - y0) / 2.0;
  const PylithScalar j3 = (z1 - z0) / 2.0;
  const PylithScalar jdet = sqrt(j1*j1 + j2*j2 + j3*j3);

  for (int i=0, iJ=0; i < npts; ++i) {
    jacobian[iJ++] = j1;
    jacobian[iJ++] = j2;
    jacobian[iJ++] = j3;
    det[i] = jdet;
  } // for

  PetscLogFlops(12);
} // jacobian


// ----------------------------------------------------------------------
// Compute minimum width across cell.
PylithScalar
pylith::feassemble::GeometryLine3D::minCellWidth(const PylithScalar* coordinatesCell,
						 const int numVertices,
						 const int spaceDim) const
{ // minCellWidth
  const int numCorners = 2;
  const int dim = 3;

  assert(numCorners == numVertices);
  assert(dim == spaceDim);

  const PylithScalar xA = coordinatesCell[0];
  const PylithScalar yA = coordinatesCell[1];
  const PylithScalar zA = coordinatesCell[2];
  const PylithScalar xB = coordinatesCell[3];
  const PylithScalar yB = coordinatesCell[4];
  const PylithScalar zB = coordinatesCell[5];
    
  const PylithScalar minWidth = sqrt(pow(xB-xA,2) + pow(yB-yA,2) + pow(zB-zA,2));

  PetscLogFlops(9);

  return minWidth;
} // minCellWidth


// End of file

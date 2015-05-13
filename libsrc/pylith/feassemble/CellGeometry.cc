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

#include "CellGeometry.hh" // implementation of class methods

#include "petsc.h" // USES PetscLogFlops

#include "pylith/utils/error.h" // USES std::logic_error
#include <cstring> // USES memcpy()

#include <iostream> // USES std::cerr

// ----------------------------------------------------------------------
// Default constructor.
pylith::feassemble::CellGeometry::CellGeometry(const ShapeEnum shape,
					       const int spaceDim) :
  _orientFn(0),
  _spaceDim(spaceDim),
  _shape(shape)
{ // constructor
  switch (shape)
    { // switch
    case LINE :
      _orientFn = _orient1D;
      break;
    case TRIANGLE :
    case QUADRILATERAL :
      _orientFn = _orient2D;
      break;
    case TETRAHEDRON :
    case HEXAHEDRON :
      break;
    default:
      std::cerr 
	<< "Could not find orientation function for cell with shape "
	<< shape << ".";
      assert(0);
      throw std::logic_error("Bad shape in CellGeometry.");
    } // switch
} // constructor

// ----------------------------------------------------------------------
// Default destructor.
pylith::feassemble::CellGeometry::~CellGeometry(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Copy constructor.
pylith::feassemble::CellGeometry::CellGeometry(const CellGeometry& g) :
  _orientFn(g._orientFn),
  _spaceDim(g._spaceDim),
  _shape(g._shape)
{ // copy constructor
} // copy constructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::feassemble::CellGeometry::deallocate(void)
{ // deallocate
} // deallocate
  
// ----------------------------------------------------------------------
// Get dimension of cell.
int
pylith::feassemble::CellGeometry::cellDim(void) const
{ // cellDim
  int dim = 0;
  switch (_shape)
    { // switch
    case LINE :
      dim = 1;
      break;
    case TRIANGLE :
    case QUADRILATERAL :
      dim = 2;
      break;
    case TETRAHEDRON :
    case HEXAHEDRON :
      dim = 3;
      break;
    default:
      std::cerr 
	<< "Could not find dimension of cell with shape " << _shape << ".";
      assert(0);
      throw std::logic_error("Bad shape in CellGeometry.");
    } // switch
  return dim;
} // cellDim

// ----------------------------------------------------------------------
// Get number of corners in cell.
int
pylith::feassemble::CellGeometry::numCorners(void) const
{ // numCorners
  int corners = 0;
  switch (_shape)
    { // switch
    case LINE :
      corners = 2;
      break;
    case TRIANGLE :
      corners = 3;
      break;
    case QUADRILATERAL :
      corners = 4;
      break;
    case TETRAHEDRON :
      corners = 4;
      break;
    case HEXAHEDRON :
      corners = 8;
      break;
    default:
      std::cerr 
	<< "Could not find number of corners for cell with shape "
	<< _shape << ".";
      assert(0);
      throw std::logic_error("Bad shape in CellGeometry.");
    } // switch
  return corners;
} // numCorners

// ----------------------------------------------------------------------
// Set coordinates of vertices in reference cell.
void
pylith::feassemble::CellGeometry::_setVertices(const PylithScalar* vertices,
					       const int numVertices,
					       const int dim)
{ // _setVertices
  assert(numCorners() == numVertices);
  assert(cellDim() == dim);
  const int nbytes = numVertices*dim*sizeof(PylithScalar);
  _vertices.resize(numVertices*dim);
  memcpy(&_vertices[0], vertices, nbytes);
} // _setVertices

// ----------------------------------------------------------------------
// Compute orientation of 1-D cell.
void
pylith::feassemble::CellGeometry::_orient1D(scalar_array* orientation,
					    const scalar_array& jacobian,
					    const PylithScalar jacobianDet,
					    const scalar_array& upDir)
{ // _orient1D
  const int orientSize = 4;
  assert(orientation);
  assert(orientSize == orientation->size());
  const int jacobianSize = 2;
  assert(jacobianSize == jacobian.size());

  // cellDim = 1
  // spaceDim = 2;

  const PylithScalar j1 = jacobian[0];
  const PylithScalar j2 = jacobian[1];
  (*orientation)[0] =  j1;
  (*orientation)[1] =  j2;
  (*orientation)[2] =  j2;
  (*orientation)[3] = -j1;
  PetscLogFlops(1);
} // _orient1D
		
// ----------------------------------------------------------------------
// Compute orientation of 2-D cell.
void
pylith::feassemble::CellGeometry::_orient2D(scalar_array* orientation,
					    const scalar_array& jacobian,
					    const PylithScalar jacobianDet,
					    const scalar_array& upDir)
{ // _orient2D
  const int orientSize = 9;
  const int jacobianSize = 6;
  assert(orientation);
  assert(orientSize == orientation->size());
  assert(jacobianSize == jacobian.size());
  assert(3 == upDir.size());

  // const int cellDim = 2;
  // const int spaceDim = 3;

  const PylithScalar j00 = jacobian[0];
  const PylithScalar j01 = jacobian[1];
  const PylithScalar j10 = jacobian[2];
  const PylithScalar j11 = jacobian[3];
  const PylithScalar j20 = jacobian[4];
  const PylithScalar j21 = jacobian[5];

  // Compute normal using Jacobian
  PylithScalar r0 =  j10*j21 - j20*j11;
  PylithScalar r1 = -j00*j21 + j20*j01;
  PylithScalar r2 =  j00*j11 - j10*j01;
  // Make unit vector
  PylithScalar mag = sqrt(r0*r0 + r1*r1 + r2*r2);
  assert(mag > 0.0);
  r0 /= mag;
  r1 /= mag;
  r2 /= mag;
  
  // Compute along-strike direction; cross product of "up" and normal
  PylithScalar p0 =  upDir[1]*r2 - upDir[2]*r1;
  PylithScalar p1 = -upDir[0]*r2 + upDir[2]*r0;
  PylithScalar p2 =  upDir[0]*r1 - upDir[1]*r0;
  // Make unit vector
  mag = sqrt(p0*p0 + p1*p1 + p2*p2);
  if (mag < 1.0e-6) {
    std::ostringstream msg;
    msg << "Error computing orientation of cell face. Cannot resolve tangential components into unambigious directions.\n"
	<< "Up direction ("
	<< upDir[0] << ", " << upDir[1] << ", " << upDir[2] << ") "
	<< " cannot be parallel to the face normal ("
	<< r0 << ", " << r1 << ", " << r2 << ").\n"
	<< "If the face is horizontal, adjust the up_dir parameter.";
    throw std::runtime_error(msg.str());
  } // if
  p0 /= mag;
  p1 /= mag;
  p2 /= mag;
  
  // Compute up-dip direction; cross product of normal and along-strike
  const PylithScalar q0 =  r1*p2 - r2*p1;
  const PylithScalar q1 = -r0*p2 + r2*p0;
  const PylithScalar q2 =  r0*p1 - r1*p0;
  mag = sqrt(q0*q0 + q1*q1 + q2*q2);
  assert(mag > 0.0);
  
  const PylithScalar wt = jacobianDet;
  (*orientation)[0] =  p0*wt;
  (*orientation)[1] =  p1*wt;
  (*orientation)[2] =  p2*wt;
  (*orientation)[3] =  q0*wt;
  (*orientation)[4] =  q1*wt;
  (*orientation)[5] =  q2*wt;
  (*orientation)[6] =  r0*wt;
  (*orientation)[7] =  r1*wt;
  (*orientation)[8] =  r2*wt;
  PetscLogFlops(63);
} // _orient2D


// End of file

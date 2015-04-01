// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Fault.hh" // implementation of object methods

#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::Fault::Fault(void) :
  _faultMesh(0),
  _id(0),
  _label(""),
  _edge("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor.
pylith::faults::Fault::~Fault(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void 
pylith::faults::Fault::deallocate(void)
{ // deallocate
  delete _faultMesh; _faultMesh = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Get mesh associated with fault fields.
const pylith::topology::Mesh&
pylith::faults::Fault::faultMesh(void) const
{ // faultMesh
  return *_faultMesh;
} // faultMesh


// ----------------------------------------------------------------------
// Get dimension of mesh.
int
pylith::faults::Fault::dimension(void) const
{ // dimension
  return (_faultMesh) ? _faultMesh->dimension() : 0;
} // dimension


// ----------------------------------------------------------------------
// Get number of vertices per cell for mesh.
int
pylith::faults::Fault::numCorners(void) const
{ // coneSize
  return (_faultMesh) ? _faultMesh->numCorners() : 0;
} // coneSize


// ----------------------------------------------------------------------
// Get number of vertices in mesh.
int
pylith::faults::Fault::numVertices(void) const
{ // numVertices
  return (_faultMesh) ? _faultMesh->numVertices() : 0;
} // numVertices


// ----------------------------------------------------------------------
// Get number of cells in mesh.
int
pylith::faults::Fault::numCells(void) const
{ // numCells
  return (_faultMesh) ? _faultMesh->numCells() : 0;
} // numCells


// End of file 

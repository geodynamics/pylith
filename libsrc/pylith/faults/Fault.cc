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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "Fault.hh" // implementation of object methods

#include "pylith/topology/SubMesh.hh" // USES SubMesh

// ----------------------------------------------------------------------
// Default constructor.
pylith::faults::Fault::Fault(void) :
  _id(0),
  _label(""),
  _faultMesh(0)
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
const pylith::topology::SubMesh&
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
// Get representative cone size for mesh.
int
pylith::faults::Fault::coneSize(void) const
{ // coneSize
  
  return (_faultMesh && numCells() > 0) ? 
    _faultMesh->sieveMesh()->getSieve()->getConeSize(*_faultMesh->sieveMesh()->heightStratum(1)->begin()) : 0;
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
  return (_faultMesh && !_faultMesh->sieveMesh().isNull() && _faultMesh->sieveMesh()->height() > 0) ? 
    _faultMesh->sieveMesh()->heightStratum(1)->size() : 0;
} // numCells


// End of file 

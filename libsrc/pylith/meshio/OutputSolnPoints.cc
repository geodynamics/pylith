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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "OutputSolnPoints.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnPoints::OutputSolnPoints(void) :
  _mesh(0),
  _pointsMesh(0),
  _interpolator(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnPoints::~OutputSolnPoints(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnPoints::deallocate(void)
{ // deallocate
  OutputManager<topology::Mesh, topology::Field<topology::Mesh> >::deallocate();

  if (_interpolator) {
    assert(_mesh);
    DMMeshInterpolationDestroy(_mesh->sieveMesh(), _interpolator);
    delete _interpolater; _interpolator = 0;
  } // if

  delete _mesh; _mesh = 0;
  delete _pointsMesh; _pointsMesh = 0;
} // deallocate
  
// ----------------------------------------------------------------------
// Get mesh associated with points.
const pylith::topology::Mesh&
pylith::meshio::OutputSolnPoints::createPointsMesh(const PylithScalar* points,
						   const int numPoints,
						   const int spaceDim)
{ // createPointsMesh
  //  delete _pointsMesh; _pointsMesh = new topology::Mesh();
  assert(_pointsMesh);
  return *_pointsMesh;
} // createPointsMesh


// End of file 

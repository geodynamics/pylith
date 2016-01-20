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

#include "OutputSolnSubset.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::OutputSolnSubset::OutputSolnSubset(void) :
  _label(""),
  _submesh(0)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::OutputSolnSubset::~OutputSolnSubset(void)
{ // destructor
  deallocate();
} // destructor  

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::OutputSolnSubset::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  OutputManager::deallocate();

  delete _submesh; _submesh = 0;

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set label identifier for subdomain.
void
pylith::meshio::OutputSolnSubset::label(const char* value)
{ // label
  PYLITH_METHOD_BEGIN;

  _label = value;

  PYLITH_METHOD_END;
} // label

// ----------------------------------------------------------------------
// Verify configuration is acceptable.
void
pylith::meshio::OutputSolnSubset::verifyConfiguration(const topology::Mesh& mesh) const
{ // verifyConfiguration
  PYLITH_METHOD_BEGIN;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscBool hasLabel = PETSC_FALSE;
  PetscErrorCode err = DMHasLabel(dmMesh, _label.c_str(), &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Mesh missing group of vertices '" << _label
	<< " for subdomain output.";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // verifyConfiguration

// ----------------------------------------------------------------------
// Get mesh associated with subdomain.
const pylith::topology::Mesh&
pylith::meshio::OutputSolnSubset::subdomainMesh(const topology::Mesh& mesh)
{ // subdomainMesh
  PYLITH_METHOD_BEGIN;

  delete _submesh; _submesh = new topology::Mesh(mesh, _label.c_str());assert(_submesh);

  PYLITH_METHOD_RETURN(*_submesh);
} // subdomainMesh


// End of file 

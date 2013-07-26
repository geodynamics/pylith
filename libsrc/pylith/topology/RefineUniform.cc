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
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

#include <portinfo>

#include "RefineUniform.hh" // implementation of class methods

#include "Mesh.hh" // USES Mesh

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ----------------------------------------------------------------------
// Constructor
pylith::topology::RefineUniform::RefineUniform(void)
{ // constructor
} // constructor
 
// ----------------------------------------------------------------------
// Destructor
pylith::topology::RefineUniform::~RefineUniform(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::topology::RefineUniform::deallocate(void)
{ // deallocate
} // deallocate

// ----------------------------------------------------------------------
// Refine mesh.
void
pylith::topology::RefineUniform::refine(Mesh* const newMesh,
					const Mesh& mesh,
					const int levels)
{ // refine
  assert(newMesh);

  PetscErrorCode err;

  PetscDM dmNew = newMesh->dmMesh();
  PetscDM dmOrig = mesh.dmMesh();
  
  PetscInt meshDepth = 0;
  err = DMPlexGetDepth(dmOrig, &meshDepth);

  const int meshDim = mesh.dimension();
  if (meshDim > 0 && meshDepth !=  meshDim) {
    std::ostringstream msg;
    msg << "Mesh refinement for uninterpolated meshes not supported.\n"
	<< "Turn on interpolated meshes using 'interpolate' mesh generator property.";
    throw std::runtime_error(msg.str());
  } // if

  err = DMPlexSetRefinementUniform(dmOrig, PETSC_TRUE);PYLITH_CHECK_ERROR(err);

  for (int i=0; i < levels; ++i) {
    err = DMRefine(dmOrig, mesh.comm(), &dmNew);PYLITH_CHECK_ERROR(err);
  } // for

  // newMesh->view("REFINED MESH");
} // refine
    

// End of file 

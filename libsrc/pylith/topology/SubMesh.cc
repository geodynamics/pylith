
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

#include "SubMesh.hh" // implementation of class methods

#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional
#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/utils/petscfwd.h" // USES PetscVec
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include <Selection.hh> // USES ALE::Selection

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::SubMesh::SubMesh(void) :
  _newMesh(NULL),
  _coordsys(0),
  _debug(false)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with mesh and label for vertices marking boundary.
pylith::topology::SubMesh::SubMesh(const Mesh& mesh,
				   const char* label) :
  _newMesh(NULL),
  _coordsys(0),
  _debug(false)
{ // constructor
  PYLITH_METHOD_BEGIN;

  createSubMesh(mesh, label);

  PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::SubMesh::~SubMesh(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::SubMesh::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _coordsys; _coordsys = 0;
  PetscErrorCode err = DMDestroy(&_newMesh);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Create PETSc mesh.
void
pylith::topology::SubMesh::createSubMesh(const Mesh& mesh,
					 const char* label)
{ // createSubMesh
  PYLITH_METHOD_BEGIN;

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscBool hasLabel;
  PetscErrorCode err;

  err = DMPlexHasLabel(dmMesh, label, &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in DM mesh.";
    throw std::runtime_error(msg.str());
  } // if

  /* TODO: Add creation of pointSF for submesh */
  err = DMDestroy(&_newMesh);PYLITH_CHECK_ERROR(err);
  err = DMPlexCreateSubmesh(dmMesh, label, 1, &_newMesh);PYLITH_CHECK_ERROR(err);

  // Set data from mesh.
  coordsys(mesh);

  // Set name
  std::string meshLabel = "subdomain_" + std::string(label);
  err = PetscObjectSetName((PetscObject) _newMesh, meshLabel.c_str());PYLITH_CHECK_ERROR(err);
  PetscInt maxConeSizeLocal, maxConeSize = 0;

  err = DMPlexGetMaxSizes(dmMesh, &maxConeSizeLocal, NULL);PYLITH_CHECK_ERROR(err);
  err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                      PetscObjectComm((PetscObject) dmMesh)); PYLITH_CHECK_ERROR(err);

  if (maxConeSize <= 0) {
    std::ostringstream msg;
    msg << "Error while creating submesh. Submesh '" 
	<< label << "' does not contain any cells.\n"
	<< "Submeshes must be one dimension lower than the domain mesh.";
    throw std::runtime_error(msg.str());
  } // if

  PYLITH_METHOD_END;
} // createSubMesh

// ----------------------------------------------------------------------
// Set coordinate system using mesh.
void
pylith::topology::SubMesh::coordsys(const Mesh& mesh)
{ // coordsys
  PYLITH_METHOD_BEGIN;

  delete _coordsys; _coordsys = 0;
  const spatialdata::geocoords::CoordSys* cs = mesh.coordsys();
  if (cs) {
    _coordsys = cs->clone();assert(_coordsys);
    _coordsys->initialize();
  } // if

  PYLITH_METHOD_END;
} // coordsys

// ----------------------------------------------------------------------
// Initialize the finite-element mesh.
void 
pylith::topology::SubMesh::initialize(void)
{ // initialize
  PYLITH_METHOD_BEGIN;

  if (_coordsys)
    _coordsys->initialize();

  PYLITH_METHOD_END;
} // initialize

// End of file 

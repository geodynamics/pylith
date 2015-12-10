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
#include <stdexcept>

#include "Mesh.hh" // implementation of class methods

#include "MeshOps.hh" // USES MeshOps

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/petscfwd.h" // USES PetscVec
#include "pylith/utils/error.h" // USES PYLITH_CHECK_ERROR

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(const bool isSubMesh) :
  _dmMesh(NULL),
  _numNormalCells(0),
  _numCohesiveCells(0),
  _numNormalVertices(0),
  _numShadowVertices(0),
  _numLagrangeVertices(0),
  _coordsys(0),
  _debug(false),
  _isSubMesh(isSubMesh)
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Constructor with dimension and communicator.
pylith::topology::Mesh::Mesh(const int dim,
			     const MPI_Comm& comm) :
  _dmMesh(NULL),
  _numNormalCells(0),
  _numCohesiveCells(0),
  _numNormalVertices(0),
  _numShadowVertices(0),
  _numLagrangeVertices(0),
  _coordsys(0),
  _debug(false),
  _isSubMesh(false)
{ // constructor
  PYLITH_METHOD_BEGIN;

  PetscErrorCode err;
  err = DMCreate(comm, &_dmMesh);PYLITH_CHECK_ERROR(err);
  err = DMSetType(_dmMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
  err = DMSetDimension(_dmMesh, dim);PYLITH_CHECK_ERROR(err);
  err = PetscObjectSetName((PetscObject) _dmMesh, "domain");PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // constructor

// ----------------------------------------------------------------------
// Create submesh.
pylith::topology::Mesh::Mesh(const Mesh& mesh,
			     const char* label) :
  _dmMesh(NULL),
  _numNormalCells(0),
  _numCohesiveCells(0),
  _numNormalVertices(0),
  _numShadowVertices(0),
  _numLagrangeVertices(0),
  _coordsys(0),
  _debug(mesh._debug),
  _isSubMesh(true)
{ // Submesh constructor
  PYLITH_METHOD_BEGIN;

  assert(label);

  coordsys(mesh.coordsys());

  PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
  PetscErrorCode err;

  PetscBool hasLabel = PETSC_FALSE;
  err = DMPlexHasLabel(dmMesh, label, &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Could not find group of points '" << label << "' in PETSc DM mesh.";
    throw std::runtime_error(msg.str());
  } // if

  if (mesh.dimension() < 1) {
    throw std::logic_error("INTERNAL ERROR in MeshOps::createSubMesh()\n"
			   "Cannot create submesh for mesh with dimension < 1.");
  } // if

  /* TODO: Add creation of pointSF for submesh */
  DMLabel l;
  err = DMPlexGetLabel(dmMesh, label, &l);PYLITH_CHECK_ERROR(err);
  err = DMPlexCreateSubmesh(dmMesh, l, 1, &_dmMesh);PYLITH_CHECK_ERROR(err);

  PetscInt maxConeSizeLocal = 0, maxConeSize = 0;
  err = DMPlexGetMaxSizes(_dmMesh, &maxConeSizeLocal, NULL);PYLITH_CHECK_ERROR(err);
  err = MPI_Allreduce(&maxConeSizeLocal, &maxConeSize, 1, MPI_INT, MPI_MAX,
                      PetscObjectComm((PetscObject) _dmMesh)); PYLITH_CHECK_ERROR(err);

  if (maxConeSize <= 0) {
    std::ostringstream msg;
    msg << "Error while creating submesh. Submesh '" << label << "' does not contain any cells.\n";
    throw std::runtime_error(msg.str());
  } // if

  // Set name
  std::string meshLabel = "subdomain_" + std::string(label);
  err = PetscObjectSetName((PetscObject) _dmMesh, meshLabel.c_str());PYLITH_CHECK_ERROR(err);

  // Set lengthscale
  PylithScalar lengthScale;
  err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
  err = DMPlexSetScale(_dmMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

  // Check topology
  MeshOps::checkTopology(*this);

  PYLITH_METHOD_END;
} // SubMesh constructor
		     

// ----------------------------------------------------------------------
// Default destructor
pylith::topology::Mesh::~Mesh(void)
{ // destructor
  deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Mesh::deallocate(void)
{ // deallocate
  PYLITH_METHOD_BEGIN;

  delete _coordsys; _coordsys = 0;
  PetscErrorCode err = DMDestroy(&_dmMesh);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_END;
} // deallocate
  
// ----------------------------------------------------------------------
// Set coordinate system.
void
pylith::topology::Mesh::coordsys(const spatialdata::geocoords::CoordSys* cs)
{ // coordsys
  PYLITH_METHOD_BEGIN;

  delete _coordsys; _coordsys = (cs) ? cs->clone() : 0;
  if (_coordsys)
    _coordsys->initialize();

  PYLITH_METHOD_END;
} // coordsys

// ----------------------------------------------------------------------
// Get MPI communicator associated with mesh.
MPI_Comm
pylith::topology::Mesh::comm(void) const
{ // comm
  PYLITH_METHOD_BEGIN;
  
  MPI_Comm comm;
  if (_dmMesh) {
    PetscObjectGetComm((PetscObject) _dmMesh, &comm);
  } else {
    comm = PETSC_COMM_WORLD;
  } // if/else
  
  PYLITH_METHOD_RETURN(comm);
} // comm
    
// ----------------------------------------------------------------------
// Get MPI rank.
int
pylith::topology::Mesh::commRank(void) const
{ // commRank
  PYLITH_METHOD_BEGIN;

  int rank = 0;
  MPI_Comm comm;
  if (_dmMesh) {
    PetscObjectGetComm((PetscObject) _dmMesh, &comm);
  } else {
    comm = PETSC_COMM_WORLD;
  } // if/else
  MPI_Comm_rank(comm, &rank);

  PYLITH_METHOD_RETURN(rank);
} // commRank

// ----------------------------------------------------------------------
// Print mesh to stdout.
void
pylith::topology::Mesh::view(const char* viewOption) const
{ // view
  PYLITH_METHOD_BEGIN;

  assert(_dmMesh);

  PetscErrorCode err;
  if (strlen(viewOption) > 0) {
    const char* label = 0;
    err = PetscObjectGetName((PetscObject) _dmMesh, &label);PYLITH_CHECK_ERROR(err);

    std::ostringstream optionname, optionprefix;
    optionprefix << label << "_";
    optionname  << "-" << label << "_dm_view";

    err = DMSetOptionsPrefix(_dmMesh, optionprefix.str().c_str());PYLITH_CHECK_ERROR(err);
    err = PetscOptionsSetValue(NULL, optionname.str().c_str(), viewOption);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(_dmMesh, NULL, "-dm_view");PYLITH_CHECK_ERROR(err);

  } else {
    err = DMView(_dmMesh, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
  } // if/else

  PYLITH_METHOD_END;
} // view

// ----------------------------------------------------------------------
// Return the names of all vertex groups.
void
pylith::topology::Mesh::groups(int* numNames, 
			       char*** names) const
{ // groups
  PYLITH_METHOD_BEGIN;

  assert(numNames);
  assert(names);
  *numNames = 0;
  *names = 0;

  if (_dmMesh) {
    PetscErrorCode err = 0;

    PetscInt numLabels = 0;
    err = DMPlexGetNumLabels(_dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);

    *numNames = numLabels;
    *names = new char*[numLabels];
    for (int iLabel=0; iLabel < numLabels; ++iLabel) {
      const char* namestr = NULL;
      err = DMPlexGetLabelName(_dmMesh, iLabel, &namestr);PYLITH_CHECK_ERROR(err);
      // Must return char* that SWIG can deallocate.
      const char len = strlen(namestr);
      char* newName = 0;
      if (len > 0) {
	newName = new char[len+1];
	strncpy(newName, namestr, len+1);
      } else {
	newName = new char[1];
	newName[0] ='\0';
      } // if/else
      (*names)[iLabel] = newName;
    } // for
  } // if

  PYLITH_METHOD_END;
} // groups

// ----------------------------------------------------------------------
// Return the names of all vertex groups.
int
pylith::topology::Mesh::groupSize(const char *name)
{ // groupSize
  PYLITH_METHOD_BEGIN;

  assert(_dmMesh);

  PetscErrorCode err = 0;

  PetscBool hasLabel = PETSC_FALSE;
  err = DMPlexHasLabel(_dmMesh, name, &hasLabel);PYLITH_CHECK_ERROR(err);
  if (!hasLabel) {
    std::ostringstream msg;
    msg << "Cannot get size of group '" << name << "'. Group missing from mesh.";
    throw std::runtime_error(msg.str());
  } // if

  PetscInt size = 0;
  err = DMPlexGetLabelSize(_dmMesh, name, &size);PYLITH_CHECK_ERROR(err);

  PYLITH_METHOD_RETURN(size);
} // groupSize


// End of file 

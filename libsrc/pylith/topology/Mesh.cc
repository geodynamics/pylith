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
// Copyright (c) 2010-2017 University of California, Davis
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
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(const bool isSubmesh) :
    _dmMesh(NULL),
    _coordsys(0),
    _debug(false),
    _isSubmesh(isSubmesh),
    _isSimplex(true) {}


// ----------------------------------------------------------------------
// Constructor with dimension and communicator.
pylith::topology::Mesh::Mesh(const int dim,
                             const MPI_Comm& comm) :
    _dmMesh(NULL),
    _coordsys(0),
    _debug(false),
    _isSubmesh(false),
    _isSimplex(true) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    err = DMCreate(comm, &_dmMesh);PYLITH_CHECK_ERROR(err);
    err = DMSetType(_dmMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
    err = DMSetDimension(_dmMesh, dim);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _dmMesh, "domain");PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // constructor


// ----------------------------------------------------------------------
// Default destructor
pylith::topology::Mesh::~Mesh(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Mesh::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _coordsys;_coordsys = 0;
    PetscErrorCode err = DMDestroy(&_dmMesh);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ----------------------------------------------------------------------
// Set DMPlex mesh.
void
pylith::topology::Mesh::dmMesh(PetscDM dm,
                               const char* label) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    err = DMDestroy(&_dmMesh);PYLITH_CHECK_ERROR(err);
    _dmMesh = dm;
    err = PetscObjectSetName((PetscObject) _dmMesh, label);PYLITH_CHECK_ERROR(err);

    _isSimplex = MeshOps::isSimplexMesh(*this);

    PYLITH_METHOD_END;
}


// ----------------------------------------------------------------------
// Set coordinate system.
void
pylith::topology::Mesh::coordsys(const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;

    delete _coordsys;_coordsys = (cs) ? cs->clone() : 0;
    if (_coordsys) {
        _coordsys->initialize();
    }

    PYLITH_METHOD_END;
} // coordsys


// ----------------------------------------------------------------------
// Get MPI communicator associated with mesh.
MPI_Comm
pylith::topology::Mesh::comm(void) const {
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
pylith::topology::Mesh::commRank(void) const {
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
pylith::topology::Mesh::view(const char* viewOption) const {
    PYLITH_METHOD_BEGIN;

    assert(_dmMesh);

    PetscErrorCode err;
    if (strlen(viewOption) > 0) {
        const char* label = 0;
        err = PetscObjectGetName((PetscObject) _dmMesh, &label);PYLITH_CHECK_ERROR(err);

        std::ostringstream optionname, optionprefix;
        optionprefix << label << "_";
        optionname << "-" << label << "_dm_view";

        err = DMSetOptionsPrefix(_dmMesh, optionprefix.str().c_str());PYLITH_CHECK_ERROR(err);
        err = PetscOptionsSetValue(NULL, optionname.str().c_str(), viewOption);PYLITH_CHECK_ERROR(err);
        err = DMViewFromOptions(_dmMesh, NULL, "-dm_view");PYLITH_CHECK_ERROR(err);

    } else {
        err = DMView(_dmMesh, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // view


// End of file

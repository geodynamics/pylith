// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/Mesh.hh" // implementation of class methods

#include "pylith/topology/MeshOps.hh" // USES MeshOps

#include "spatialdata/geocoords/CoordSys.hh" // USES CoordSys
#include "pylith/utils/array.hh" // USES scalar_array
#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR
#include "pylith/utils/petscfwd.h" // USES PetscVec

#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
const char* pylith::topology::Mesh::cells_label_name = "material-id";

// ------------------------------------------------------------------------------------------------
// Default constructor
pylith::topology::Mesh::Mesh(void) :
    _coordSys(NULL),
    _dm(NULL) {}


// ------------------------------------------------------------------------------------------------
// Constructor with dimension and communicator.
pylith::topology::Mesh::Mesh(const int dim,
                             const MPI_Comm& comm) :
    _coordSys(NULL),
    _dm(NULL) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err;
    err = DMCreate(comm, &_dm);PYLITH_CHECK_ERROR(err);
    err = DMSetType(_dm, DMPLEX);PYLITH_CHECK_ERROR(err);
    err = DMSetDimension(_dm, dim);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) _dm, "domain");PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // constructor


// ------------------------------------------------------------------------------------------------
// Default destructor
pylith::topology::Mesh::~Mesh(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::topology::Mesh::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    delete _coordSys;_coordSys = NULL;
    PetscErrorCode err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Create clone.
pylith::topology::Mesh*
pylith::topology::Mesh::clone(void) const {
    PYLITH_METHOD_BEGIN;

    Mesh* mesh = new Mesh();assert(mesh);
    mesh->setCoordSys(this->getCoordSys());

    PetscErrorCode err;
    if (this->_dm) {
        err = DMClone(this->_dm, &mesh->_dm);

        const char* name = NULL;
        err = PetscObjectGetName((PetscObject)this->_dm, &name);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject)mesh->_dm,  name);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(mesh);
}


// ----------------------------------------------------------------------
// Get DMPlex mesh.
PetscDM
pylith::topology::Mesh::getDM(void) const {
    return _dm;
}


// ------------------------------------------------------------------------------------------------
// Set DMPlex mesh.
void
pylith::topology::Mesh::setDM(PetscDM dm,
                              const char* label) {
    PYLITH_METHOD_BEGIN;

    PetscErrorCode err = PETSC_SUCCESS;
    err = DMDestroy(&_dm);PYLITH_CHECK_ERROR(err);
    _dm = dm;
    if (label) {
        err = PetscObjectSetName((PetscObject) _dm, label);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_END;
}


// ----------------------------------------------------------------------
// Get coordinate system.
const spatialdata::geocoords::CoordSys*
pylith::topology::Mesh::getCoordSys(void) const {
    return _coordSys;
}


// ------------------------------------------------------------------------------------------------
// Set coordinate system.
void
pylith::topology::Mesh::setCoordSys(const spatialdata::geocoords::CoordSys* cs) {
    PYLITH_METHOD_BEGIN;

    delete _coordSys;_coordSys = (cs) ? cs->clone() : NULL;

    PYLITH_METHOD_END;
} // setCoordSys


// ----------------------------------------------------------------------
// Get dimension of mesh.
int
pylith::topology::Mesh::getDimension(void) const {
    PYLITH_METHOD_BEGIN;

    PetscInt dim = 0;
    if (_dm) {
        PetscErrorCode err = DMGetDimension(_dm, &dim);PYLITH_CHECK_ERROR(err);
    } // if

    PYLITH_METHOD_RETURN(dim);
}


// ------------------------------------------------------------------------------------------------
// Get MPI communicator associated with mesh.
MPI_Comm
pylith::topology::Mesh::getComm(void) const {
    PYLITH_METHOD_BEGIN;

    MPI_Comm comm;
    if (_dm) {
        PetscObjectGetComm((PetscObject) _dm, &comm);
    } else {
        comm = PETSC_COMM_WORLD;
    } // if/else

    PYLITH_METHOD_RETURN(comm);
} // comm


// ------------------------------------------------------------------------------------------------
// Get MPI rank.
int
pylith::topology::Mesh::getCommRank(void) const {
    PYLITH_METHOD_BEGIN;

    int rank = 0;
    MPI_Comm comm = this->getComm();
    MPI_Comm_rank(comm, &rank);

    PYLITH_METHOD_RETURN(rank);
} // commRank


// ------------------------------------------------------------------------------------------------
// Print mesh to stdout.
void
pylith::topology::Mesh::view(const char* viewOption) const {
    PYLITH_METHOD_BEGIN;

    assert(_dm);

    PetscErrorCode err;
    if (strlen(viewOption) > 0) {
        const char* label = 0;
        err = PetscObjectGetName((PetscObject) _dm, &label);PYLITH_CHECK_ERROR(err);

        std::ostringstream optionname, optionprefix;
        optionprefix << label << "_";
        optionname << "-" << label << "_dm_view";

        err = DMSetOptionsPrefix(_dm, optionprefix.str().c_str());PYLITH_CHECK_ERROR(err);
        err = PetscOptionsSetValue(NULL, optionname.str().c_str(), viewOption);PYLITH_CHECK_ERROR(err);
        err = DMViewFromOptions(_dm, NULL, "-dm_view");PYLITH_CHECK_ERROR(err);

    } else {
        err = DMView(_dm, PETSC_VIEWER_STDOUT_WORLD);PYLITH_CHECK_ERROR(err);
    } // if/else

    PYLITH_METHOD_END;
} // view


// End of file

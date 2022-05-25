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

#include "MeshIOPetsc.hh" // implementation of class methods

#include "MeshBuilder.hh" // USES MeshBuilder
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "spatialdata/utils/LineParser.hh" // USES LineParser

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*

#include <iomanip> // USES setw(), setiosflags(), resetiosflags()
#include <strings.h> // USES strcasecmp()
#include <cassert> // USES assert()
#include <fstream> // USES std::ifstream, std::ofstream
#include <stdexcept> // USES std::runtime_error
#include <sstream> // USES std::ostringstream
#include <typeinfo> // USES std::typeid

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace meshio {
        class _MeshIOPetsc {
public:

            /** Remove everything but cells from label for materials.
             *
             * @param[inout] dmMesh PETSc DM for mesh.
             */
            static
            void fixMaterialLabel(PetscDM* dmMesh);

        }; // _MeshIOPetsc
    } // meshio
} // pylith

// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOPetsc::MeshIOPetsc(void) :
    _filename(""),
    _prefix("") {
    PyreComponent::setName("meshiopetsc");
} // constructor


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::meshio::MeshIOPetsc::~MeshIOPetsc(void) {
    deallocate();
} // destructor


// ------------------------------------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::MeshIOPetsc::deallocate(void) {
    PYLITH_METHOD_BEGIN;

    MeshIO::deallocate();

    PYLITH_METHOD_END;
} // deallocate


// ------------------------------------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOPetsc::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");
    assert(_mesh);

    const size_t noptions = 3;
    std::string options[noptions*2] = {
        "-" + _prefix + "dm_plex_filename", _filename,
        "-" + _prefix + "dm_plex_gmsh_use_regions", "",
        "-" + _prefix + "dm_plex_gmsh_mark_vertices", "",
    };

    PetscErrorCode err;
    if (!_filename.empty()) {
        for (size_t i = 0; i < noptions; ++i) {
            err = PetscOptionsSetValue(NULL, options[2*i+0].c_str(), options[2*i+1].c_str());
        } // for
    } // if

    PetscDM dmMesh = NULL;
    err = DMCreate(_mesh->getComm(), &dmMesh);PYLITH_CHECK_ERROR(err);
    err = DMSetType(dmMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
    if (!_prefix.empty()) {
        err = PetscObjectSetOptionsPrefix((PetscObject) dmMesh, _prefix.c_str());PYLITH_CHECK_ERROR(err);
    } // if
    err = DMPlexDistributeSetDefault(dmMesh, PETSC_FALSE);PYLITH_CHECK_ERROR(err);
    err = DMSetFromOptions(dmMesh);PYLITH_CHECK_ERROR(err);
    _MeshIOPetsc::fixMaterialLabel(&dmMesh);
    _mesh->setDM(dmMesh);

    PYLITH_METHOD_END;
} // read


// ------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOPetsc::_write(void) const {  }


// ------------------------------------------------------------------------------------------------
// Fix material label;
void
pylith::meshio::_MeshIOPetsc::fixMaterialLabel(PetscDM* dmMesh) {
    PYLITH_METHOD_BEGIN;
    assert(dmMesh);
    PetscErrorCode err = 0;
    const char* const labelName = pylith::topology::Mesh::cells_label_name;

    const PetscInt cellHeight = 0;
    PetscInt cStart = -1, cEnd = -1;
    err = DMPlexGetHeightStratum(*dmMesh, cellHeight, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    if (cStart == cEnd) {
        err = DMCreateLabel(*dmMesh, labelName);PYLITH_CHECK_ERROR(err);
        PYLITH_METHOD_END;
    } // if

    PetscDMLabel dmLabel = NULL;
    PetscInt pStart = -1, pEnd = -1;
    err = DMGetLabel(*dmMesh, labelName, &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetBounds(dmLabel, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    if (pStart == cStart) { pStart = cEnd; }
    if (pEnd == cEnd) { pEnd = cStart; }

    for (PetscInt point = pStart; point < pEnd; ++point) {
        PetscBool hasLabel = PETSC_FALSE;
        err = DMLabelHasPoint(dmLabel, point, &hasLabel);PYLITH_CHECK_ERROR(err);
        if (hasLabel) {
            PetscInt labelValue;
            err = DMLabelGetValue(dmLabel, point, &labelValue);PYLITH_CHECK_ERROR(err);
            err = DMLabelClearValue(dmLabel, point, labelValue);PYLITH_CHECK_ERROR(err);
        } // if
    } // for

    PYLITH_METHOD_END;
} // fixMaterialLabel

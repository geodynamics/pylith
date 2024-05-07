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

#include "pylith/meshio/MeshIOPetsc.hh" // implementation of class methods

#include "pylith/meshio/MeshBuilder.hh" // USES MeshBuilder
#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // USES scalar_array, int_array, string_vector
#include "spatialdata/utils/LineParser.hh" // USES LineParser

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/journals.hh" // USES PYLITH_COMPONENT_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include "petscviewerhdf5.h"
#include "petscdmplex.h"

#include <set> // USES std::set
#include <cassert> // USES assert()
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

            /** Create labels for boundaries based on vertices.
             *
             * For each boundary condition label:
             * 1. Remove all points that are not vertices from label.
             * 2. Add edges and faces to label based upon vertices.
             *
             * @param[inout] dmMesh PETSc DM for mesh.
             */
            static
            void fixBoundaryLabels(PetscDM* dmMesh);

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt read;
                static PylithInt fixMaterialLabel;
                static PylithInt fixBoundaryLabels;
            };

        }; // _MeshIOPetsc
    } // meshio
} // pylith

pylith::utils::EventLogger pylith::meshio::_MeshIOPetsc::Events::logger;
PylithInt pylith::meshio::_MeshIOPetsc::Events::read;
PylithInt pylith::meshio::_MeshIOPetsc::Events::fixMaterialLabel;
PylithInt pylith::meshio::_MeshIOPetsc::Events::fixBoundaryLabels;

// ------------------------------------------------------------------------------------------------
void
pylith::meshio::_MeshIOPetsc::Events::init(void) {
    logger.setClassName("MeshIOPetsc");
    logger.initialize();
    read = logger.registerEvent("PL:MeshIOPetsc:read");
    fixMaterialLabel = logger.registerEvent("PL:MeshIOPetsc:fixMaterialLabel");
    fixBoundaryLabels = logger.registerEvent("PL:MeshIOPetsc:fixBoundaryLabels");
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOPetsc::MeshIOPetsc(void) :
    _filename(""),
    _prefix(""),
    _format(HDF5) {
    PyreComponent::setName("meshiopetsc");
    _MeshIOPetsc::Events::init();
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
// Set filename for PETSc.
void
pylith::meshio::MeshIOPetsc::setFilename(const char* name) {
    _filename = name;

    const std::string gmshSuffix = ".msh";
    const std::string hdf5Suffix = ".h5";
    if (gmshSuffix == _filename.substr(_filename.size()-gmshSuffix.size(), gmshSuffix.size())) {
        _format = GMSH;
    } else if (hdf5Suffix == _filename.substr(_filename.size()-hdf5Suffix.size(), hdf5Suffix.size())) {
        _format = HDF5;
    } else {
        PYLITH_COMPONENT_LOGICERROR("Could not determine format for mesh file " << _filename << ".");
    } // if/else
}


// ------------------------------------------------------------------------------------------------
// Get filename for PETSc.
const char*
pylith::meshio::MeshIOPetsc::getFilename(void) const {
    return _filename.c_str();
}


// ------------------------------------------------------------------------------------------------
// Set options prefix for mesh.
void
pylith::meshio::MeshIOPetsc::setPrefix(const char* name) {
    _prefix = name;
}


// ------------------------------------------------------------------------------------------------
// Get options prefix for mesh.
const char*
pylith::meshio::MeshIOPetsc::getPrefix(void) const {
    return _prefix.c_str();
}


// ------------------------------------------------------------------------------------------------
// Set mesh format.
void
pylith::meshio::MeshIOPetsc::setFormat(Format value) {
    _format = value;
}


// ------------------------------------------------------------------------------------------------
// Get mesh format.
pylith::meshio::MeshIOPetsc::Format
pylith::meshio::MeshIOPetsc::getFormat(void) const {
    return _format;
}


// ------------------------------------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOPetsc::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");
    _MeshIOPetsc::Events::logger.eventBegin(_MeshIOPetsc::Events::read);
    assert(_mesh);

    PetscErrorCode err = PETSC_SUCCESS;
    if (!_filename.empty()) {
        size_t noptions = 1;
        pylith::string_vector options(noptions*2);
        options[0] = "-" + _prefix + "dm_plex_filename";
        options[1] = _filename;
        if (GMSH == _format) {
            noptions += 2;
            options.resize(noptions*2);
            options[2] = "-" + _prefix + "dm_plex_gmsh_use_regions";
            options[3] = "";

            options[4] = "-" + _prefix + "dm_plex_gmsh_mark_vertices";
            options[5] = "";
        } // if

        for (size_t i = 0; i < noptions; ++i) {
            err = PetscOptionsSetValue(NULL, options[2*i+0].c_str(), options[2*i+1].c_str());
        } // for
    } // if

    PetscDM dmMesh = NULL;
    try {
        err = DMCreate(_mesh->getComm(), &dmMesh);PYLITH_CHECK_ERROR(err);
        err = DMSetType(dmMesh, DMPLEX);PYLITH_CHECK_ERROR(err);
        err = PetscObjectSetName((PetscObject) dmMesh, "domain");
        if (!_prefix.empty()) {
            err = PetscObjectSetOptionsPrefix((PetscObject) dmMesh, _prefix.c_str());PYLITH_CHECK_ERROR(err);
        } // if
        err = DMPlexDistributeSetDefault(dmMesh, PETSC_FALSE);PYLITH_CHECK_ERROR(err);
        err = DMSetFromOptions(dmMesh);PYLITH_CHECK_ERROR_MSG(err, "Error creating mesh with MeshIOPetsc.");

        _MeshIOPetsc::fixMaterialLabel(&dmMesh);
        _MeshIOPetsc::fixBoundaryLabels(&dmMesh);
        _mesh->setDM(dmMesh, "domain");
    } catch (...) {
        DMDestroy(&dmMesh);
        throw;
    } // try/catch

    _MeshIOPetsc::Events::logger.eventEnd(_MeshIOPetsc::Events::read);
    PYLITH_METHOD_END;
} // read


// ------------------------------------------------------------------------------------------------
// Write mesh to file.
void
pylith::meshio::MeshIOPetsc::_write(void) const {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_write()");
    assert(_mesh);

    PetscErrorCode err = PETSC_SUCCESS;
    PetscViewer viewer = PETSC_NULLPTR;
    if (_format == HDF5) {
        err = PetscViewerHDF5Open(PETSC_COMM_WORLD, _filename.c_str(), FILE_MODE_WRITE, &viewer);

        DMPlexStorageVersion storageVersion = PETSC_NULLPTR;
        err = PetscNew(&storageVersion);PYLITH_CHECK_ERROR(err);assert(storageVersion);
        storageVersion->major = 3;
        storageVersion->minor = 0;
        storageVersion->subminor = 0;
        err = PetscViewerHDF5SetDMPlexStorageVersionWriting(viewer, storageVersion);PYLITH_CHECK_ERROR(err);
        err = PetscFree(storageVersion);PYLITH_CHECK_ERROR(err);
    } else {
        PYLITH_COMPONENT_LOGICERROR("Unknown mesh format " << _format << ".");
    } // if/else
    err = DMView(_mesh->getDM(), viewer);PYLITH_CHECK_ERROR(err);
    err = PetscViewerDestroy(&viewer);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // _write


// ------------------------------------------------------------------------------------------------
// Remove everything but cells from label for materials.
void
pylith::meshio::_MeshIOPetsc::fixMaterialLabel(PetscDM* dmMesh) {
    PYLITH_METHOD_BEGIN;
    _MeshIOPetsc::Events::logger.eventBegin(_MeshIOPetsc::Events::fixMaterialLabel);

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

    /** Get all values for a label up front and then clear all label values for
     * all non-cell points. This is much faster than checking to see if a point
     * is in the label, getting the label value, and clearing that value.
     */
    PetscIS valuesIS = PETSC_NULLPTR;
    err = DMLabelGetNonEmptyStratumValuesIS(dmLabel, &valuesIS);PYLITH_CHECK_ERROR(err);
    PetscInt numValues;
    err = ISGetLocalSize(valuesIS, &numValues);PYLITH_CHECK_ERROR(err);
    const PetscInt* valuesIndices = PETSC_NULLPTR;
    err = ISGetIndices(valuesIS, &valuesIndices);
    for (PetscInt point = pStart; point < pEnd; ++point) {
        for (PetscInt iValue = 0; iValue < numValues; ++iValue) {
            err = DMLabelClearValue(dmLabel, point, valuesIndices[iValue]);PYLITH_CHECK_ERROR(err);
        } // for
    } // for
    err = ISRestoreIndices(valuesIS, &valuesIndices);
    err = ISDestroy(&valuesIS);

    _MeshIOPetsc::Events::logger.eventEnd(_MeshIOPetsc::Events::fixMaterialLabel);
    PYLITH_METHOD_END;
} // fixMaterialLabel


// ------------------------------------------------------------------------------------------------
// Fix boundary condition labels. Remove points other than vertices from all labels other
// than the material label.
void
pylith::meshio::_MeshIOPetsc::fixBoundaryLabels(PetscDM* dmMesh) {
    PYLITH_METHOD_BEGIN;
    _MeshIOPetsc::Events::logger.eventBegin(_MeshIOPetsc::Events::fixBoundaryLabels);
    assert(dmMesh);
    PetscErrorCode err = 0;

    // Create set with labels to ignore.
    std::set<std::string> labelsIgnore;
    labelsIgnore.insert(std::string(pylith::topology::Mesh::cells_label_name));
    labelsIgnore.insert("celltype");
    labelsIgnore.insert("depth");

    PetscInt vStart = -1, vEnd = -1;
    PetscInt eStart = -1, eEnd = -1;
    PetscInt fStart = -1, fEnd = -1;
    PetscInt cStart = -1, cEnd = -1;
    err = DMPlexGetDepthStratum(*dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(*dmMesh, 1, &eStart, &eEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(*dmMesh, 2, &fStart, &fEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(*dmMesh, 3, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);

    PetscInt numLabels = 0;
    err = DMGetNumLabels(*dmMesh, &numLabels);PYLITH_CHECK_ERROR(err);
    for (PetscInt iLabel = 0; iLabel < numLabels; ++iLabel) {
        const char* labelName = NULL;
        err = DMGetLabelName(*dmMesh, iLabel, &labelName);PYLITH_CHECK_ERROR(err);
        if (labelsIgnore.count(std::string(labelName)) > 0) { continue; }

        PetscDMLabel dmLabel = NULL;
        err = DMGetLabelByNum(*dmMesh, iLabel, &dmLabel);PYLITH_CHECK_ERROR(err);
        PetscInt pStart = -1, pEnd = -1;
        err = DMLabelGetBounds(dmLabel, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);

        bool hasEdges = false;
        bool hasFaces = false;
        bool hasCells = false;

        for (PetscInt point = pStart; point < pEnd; ++point) {
            if ((point >= vStart) && (point < vEnd)) {
                continue; // always keep vertices in label
            } // if
            if ((point >= eStart) && (point < eEnd)) {
                hasEdges = true;
            } // if
            if ((point >= fStart) && (point < fEnd)) {
                hasFaces = true;
            } // if
            if ((point >= cStart) && (point < cEnd)) {
                hasCells = true;
            } // if
            PetscBool hasLabel = PETSC_FALSE;
            err = DMLabelHasPoint(dmLabel, point, &hasLabel);PYLITH_CHECK_ERROR(err);
            if (hasLabel) {
                PetscInt labelValue;
                err = DMLabelGetValue(dmLabel, point, &labelValue);PYLITH_CHECK_ERROR(err);
                err = DMLabelClearValue(dmLabel, point, labelValue);PYLITH_CHECK_ERROR(err);
            } // if
        } // for
        err = DMLabelDestroyIndex(dmLabel);PYLITH_CHECK_ERROR(err);

        // Re-add edges, faces, and cells based upon vertices if present in original label.
        err = DMLabelComputeIndex(dmLabel);PYLITH_CHECK_ERROR(err);
        err = DMLabelGetBounds(dmLabel, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
        err = DMLabelDestroyIndex(dmLabel);PYLITH_CHECK_ERROR(err);
        for (PetscInt vertex = pStart; vertex < pEnd; ++vertex) {
            PetscInt labelValue = -1;
            err = DMLabelGetValue(dmLabel, vertex, &labelValue);PYLITH_CHECK_ERROR(err);
            if (labelValue < 1) { continue; }

            PetscInt *star = NULL, starSize;
            err = DMPlexGetTransitiveClosure(*dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
            for (PetscInt s = 0; s < starSize*2; s += 2) {
                const PetscInt point = star[s];

                // Ignore depths not in original label.
                if (!hasEdges && ((point >= eStart) && (point < eEnd))) { continue; }
                if (!hasFaces && ((point >= fStart) && (point < fEnd))) { continue; }
                if (!hasCells && ((point >= cStart) && (point < cEnd))) { continue; }

                // All vertices in closure must be in label to add point to label.
                PetscInt *closure = NULL, closureSize, value;
                PetscBool mark = PETSC_TRUE;
                err = DMPlexGetTransitiveClosure(*dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
                for (PetscInt c = 0; c < closureSize*2; c += 2) {
                    if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                        err = DMLabelGetValue(dmLabel, closure[c], &value);PYLITH_CHECK_ERROR(err);
                        if (value != labelValue) {
                            mark = PETSC_FALSE;
                            break;
                        } // if
                    } // if
                } // for
                err = DMPlexRestoreTransitiveClosure(*dmMesh, point, PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
                if (mark) {
                    err = DMLabelSetValue(dmLabel, point, labelValue);PYLITH_CHECK_ERROR(err);
                } // if
            } // for
            err = DMPlexRestoreTransitiveClosure(*dmMesh, vertex, PETSC_FALSE, &starSize, &star);PYLITH_CHECK_ERROR(err);
        } // for
        err = DMLabelComputeIndex(dmLabel);PYLITH_CHECK_ERROR(err);

    } // for

    _MeshIOPetsc::Events::logger.eventEnd(_MeshIOPetsc::Events::fixBoundaryLabels);
    PYLITH_METHOD_END;
} // fixBoundaryLabels

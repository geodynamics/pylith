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
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::meshio::MeshIOPetsc::MeshIOPetsc(void) :
    _filename(""),
    _prefix(""),
    _format(HDF5),
    _gmshMarkRecursive(false) {
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
        PYLITH_COMPONENT_LOGICERROR("Could not determine format for mesh file " << _filename << " from suffix (must be '.msh' for GMSH and '.h5' for HDF5).");
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
// Set flag for marking Gmsh vertices.
void
pylith::meshio::MeshIOPetsc::setGmshMarkRecursive(const bool value) {
    _gmshMarkRecursive = value;
}


// ------------------------------------------------------------------------------------------------
// Returns true if marking Gmsh vertices, otherwise false.
bool
pylith::meshio::MeshIOPetsc::getGmshMarkRecursive(void) const {
    return _gmshMarkRecursive;
}


// ------------------------------------------------------------------------------------------------
// Read mesh.
void
pylith::meshio::MeshIOPetsc::_read(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_COMPONENT_DEBUG("_read()");
    _MeshIOPetsc::Events::logger.eventBegin(_MeshIOPetsc::Events::read);
    assert(_mesh);

    if (!_filename.empty()) {
        size_t noptions = 1;
        pylith::string_vector options(noptions*2);
        options[0] = "-" + _prefix + "dm_plex_filename";
        options[1] = _filename;
        if (GMSH == _format) {
            noptions += 3;
            options.resize(noptions*2);
            options[2] = "-" + _prefix + "dm_plex_gmsh_use_regions";
            options[3] = "true";
            options[4] = "-" + _prefix + "dm_plex_gmsh_mark_vertices_strict"; // physical group with dim==0
            options[5] = "true";
            options[6] = "-" + _prefix + "dm_plex_gmsh_mark_vertices"; // faces, edges, & vertices; :DEPRECATED:
            options[7] = (_gmshMarkRecursive) ? "true" : "false";
        } // if

        for (size_t i = 0; i < noptions; ++i) {
            PylithCallPetsc(PetscOptionsSetValue(NULL, options[2*i+0].c_str(), options[2*i+1].c_str()));
        } // for
    } // if

    PetscDM dmMesh = NULL;
    try {
        PylithCallPetsc(DMCreate(_mesh->getComm(), &dmMesh));
        PylithCallPetsc(DMSetType(dmMesh, DMPLEX));
        PylithCallPetsc(PetscObjectSetName((PetscObject) dmMesh, "domain")); // Needed for reading PETSc HDF5
        if (!_prefix.empty()) {
            PylithCallPetsc(PetscObjectSetOptionsPrefix((PetscObject) dmMesh, _prefix.c_str()));
        } // if
        PylithCallPetsc(DMPlexDistributeSetDefault(dmMesh, PETSC_FALSE));
        PylithCallPetsc(DMSetFromOptions(dmMesh));
        _MeshIOPetsc::fixMaterialLabel(&dmMesh);
        if (_gmshMarkRecursive) {
            _MeshIOPetsc::fixBoundaryLabels(&dmMesh);
        } // if
        _mesh->setDM(dmMesh);
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

    PetscViewer viewer = PETSC_NULLPTR;
    if (_format == HDF5) {
        PylithCallPetsc(PetscViewerHDF5Open(PETSC_COMM_WORLD, _filename.c_str(), FILE_MODE_WRITE, &viewer));

        DMPlexStorageVersion storageVersion = PETSC_NULLPTR;
        PylithCallPetsc(PetscNew(&storageVersion));assert(storageVersion);
        storageVersion->major = 3;
        storageVersion->minor = 1;
        storageVersion->subminor = 0;
        PylithCallPetsc(PetscViewerHDF5SetDMPlexStorageVersionWriting(viewer, storageVersion));
        PylithCallPetsc(PetscFree(storageVersion));
    } else {
        PYLITH_COMPONENT_LOGICERROR("Unknown mesh format " << _format << ".");
    } // if/else
    PylithCallPetsc(DMView(_mesh->getDM(), viewer));
    PylithCallPetsc(PetscViewerDestroy(&viewer));

    PYLITH_METHOD_END;
} // _write


// ------------------------------------------------------------------------------------------------
// Remove everything but cells from label for materials.
void
pylith::meshio::_MeshIOPetsc::fixMaterialLabel(PetscDM* dmMesh) {
    PYLITH_METHOD_BEGIN;
    _MeshIOPetsc::Events::logger.eventBegin(_MeshIOPetsc::Events::fixMaterialLabel);

    assert(dmMesh);
    const char* const labelName = pylith::topology::Mesh::cells_label_name;

    const PetscInt cellHeight = 0;
    PetscInt cStart = -1, cEnd = -1;
    PylithCallPetsc(DMPlexGetHeightStratum(*dmMesh, cellHeight, &cStart, &cEnd));
    if (cStart == cEnd) {
        PylithCallPetsc(DMCreateLabel(*dmMesh, labelName));
        PYLITH_METHOD_END;
    } // if

    PetscDMLabel dmLabel = NULL;
    PetscInt pStart = -1, pEnd = -1;
    PylithCallPetsc(DMGetLabel(*dmMesh, labelName, &dmLabel));
    PylithCallPetsc(DMLabelGetBounds(dmLabel, &pStart, &pEnd));
    if (pStart == cStart) { pStart = cEnd; }
    if (pEnd == cEnd) { pEnd = cStart; }

    /** Get all values for a label up front and then clear all label values for
     * all non-cell points. This is much faster than checking to see if a point
     * is in the label, getting the label value, and clearing that value.
     */
    PetscIS valuesIS = PETSC_NULLPTR;
    PylithCallPetsc(DMLabelGetNonEmptyStratumValuesIS(dmLabel, &valuesIS));
    PetscInt numValues;
    PylithCallPetsc(ISGetLocalSize(valuesIS, &numValues));
    const PetscInt* valuesIndices = PETSC_NULLPTR;
    PylithCallPetsc(ISGetIndices(valuesIS, &valuesIndices));
    for (PetscInt point = pStart; point < pEnd; ++point) {
        for (PetscInt iValue = 0; iValue < numValues; ++iValue) {
            PylithCallPetsc(DMLabelClearValue(dmLabel, point, valuesIndices[iValue]));
        } // for
    } // for
    PylithCallPetsc(ISRestoreIndices(valuesIS, &valuesIndices));
    PylithCallPetsc(ISDestroy(&valuesIS));

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

    // Create set with labels to ignore.
    std::set<std::string> labelsIgnore;
    labelsIgnore.insert(std::string(pylith::topology::Mesh::cells_label_name));
    labelsIgnore.insert("celltype");
    labelsIgnore.insert("depth");

    PetscInt vStart = -1, vEnd = -1;
    PetscInt eStart = -1, eEnd = -1;
    PetscInt fStart = -1, fEnd = -1;
    PetscInt cStart = -1, cEnd = -1;
    PylithCallPetsc(DMPlexGetDepthStratum(*dmMesh, 0, &vStart, &vEnd));
    PylithCallPetsc(DMPlexGetDepthStratum(*dmMesh, 1, &eStart, &eEnd));
    PylithCallPetsc(DMPlexGetDepthStratum(*dmMesh, 2, &fStart, &fEnd));
    PylithCallPetsc(DMPlexGetDepthStratum(*dmMesh, 3, &cStart, &cEnd));

    PetscInt numLabels = 0;
    PylithCallPetsc(DMGetNumLabels(*dmMesh, &numLabels));
    for (PetscInt iLabel = 0; iLabel < numLabels; ++iLabel) {
        const char* labelName = NULL;
        PylithCallPetsc(DMGetLabelName(*dmMesh, iLabel, &labelName));
        if (labelsIgnore.count(std::string(labelName)) > 0) { continue; }

        PetscDMLabel dmLabel = NULL;
        PylithCallPetsc(DMGetLabelByNum(*dmMesh, iLabel, &dmLabel));
        PetscInt pStart = -1, pEnd = -1;
        PylithCallPetsc(DMLabelGetBounds(dmLabel, &pStart, &pEnd));

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
            PylithCallPetsc(DMLabelHasPoint(dmLabel, point, &hasLabel));
            if (hasLabel) {
                PetscInt labelValue;
                PylithCallPetsc(DMLabelGetValue(dmLabel, point, &labelValue));
                PylithCallPetsc(DMLabelClearValue(dmLabel, point, labelValue));
            } // if
        } // for
        PylithCallPetsc(DMLabelDestroyIndex(dmLabel));

        // Re-add edges, faces, and cells based upon vertices if present in original label.
        PylithCallPetsc(DMLabelComputeIndex(dmLabel));
        PylithCallPetsc(DMLabelGetBounds(dmLabel, &pStart, &pEnd));
        PylithCallPetsc(DMLabelDestroyIndex(dmLabel));
        for (PetscInt vertex = pStart; vertex < pEnd; ++vertex) {
            PetscInt labelValue = -1;
            PylithCallPetsc(DMLabelGetValue(dmLabel, vertex, &labelValue));
            if (labelValue < 1) { continue; }

            PetscInt *star = NULL, starSize;
            PylithCallPetsc(DMPlexGetTransitiveClosure(*dmMesh, vertex, PETSC_FALSE, &starSize, &star));
            for (PetscInt s = 0; s < starSize*2; s += 2) {
                const PetscInt point = star[s];

                // Ignore depths not in original label.
                if (!hasEdges && ((point >= eStart) && (point < eEnd))) { continue; }
                if (!hasFaces && ((point >= fStart) && (point < fEnd))) { continue; }
                if (!hasCells && ((point >= cStart) && (point < cEnd))) { continue; }

                // All vertices in closure must be in label to add point to label.
                PetscInt *closure = NULL, closureSize, value;
                PetscBool mark = PETSC_TRUE;
                PylithCallPetsc(DMPlexGetTransitiveClosure(*dmMesh, point, PETSC_TRUE, &closureSize, &closure));
                for (PetscInt c = 0; c < closureSize*2; c += 2) {
                    if ((closure[c] >= vStart) && (closure[c] < vEnd)) {
                        PylithCallPetsc(DMLabelGetValue(dmLabel, closure[c], &value));
                        if (value != labelValue) {
                            mark = PETSC_FALSE;
                            break;
                        } // if
                    } // if
                } // for
                PylithCallPetsc(DMPlexRestoreTransitiveClosure(*dmMesh, point, PETSC_TRUE, &closureSize, &closure));
                if (mark) {
                    PylithCallPetsc(DMLabelSetValue(dmLabel, point, labelValue));
                } // if
            } // for
            PylithCallPetsc(DMPlexRestoreTransitiveClosure(*dmMesh, vertex, PETSC_FALSE, &starSize, &star));
        } // for
        PylithCallPetsc(DMLabelComputeIndex(dmLabel));

    } // for

    _MeshIOPetsc::Events::logger.eventEnd(_MeshIOPetsc::Events::fixBoundaryLabels);
    PYLITH_METHOD_END;
} // fixBoundaryLabels

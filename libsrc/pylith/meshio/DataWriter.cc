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

#include "pylith/meshio/DataWriter.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES isCohesiveCell()

#include "pylith/utils/error.hh" \
    // USES PYLITH_METHOD_BEGIN/END

// Constructor
pylith::meshio::DataWriter::DataWriter(void) :
    _timeScale(1.0),
    _context(""),
    _isInfo(false),
    _isOpen(false) {}


// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriter::~DataWriter(void) {
    deallocate();
} // destructor


// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriter::deallocate(void) {
}


// ----------------------------------------------------------------------
// Set time scale for simulation time.
void
pylith::meshio::DataWriter::setTimeScale(const PylithScalar value) {
    PYLITH_METHOD_BEGIN;

    if (value <= 0.0) {
        std::ostringstream msg;
        msg << "Time scale for simulation time (" << value << " must be positive.";
        throw std::runtime_error(msg.str());
    } // if

    _timeScale = value;

    PYLITH_METHOD_END;
} // timeScale


// ----------------------------------------------------------------------
// Is data writer open, i.e., ready for openTimeStep()/closeTimeStep()?
bool
pylith::meshio::DataWriter::isOpen(void) const {
    return _isOpen;
} // isOpen


// ----------------------------------------------------------------------
// Prepare for writing files.
void
pylith::meshio::DataWriter::open(const pylith::topology::Mesh& mesh,
                                 const bool isInfo) {
    PYLITH_METHOD_BEGIN;

    _isInfo = isInfo;
    _isOpen = true;

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    const char* meshName = NULL;
    PetscObjectGetName((PetscObject) dmMesh, &meshName);

    std::ostringstream s;
    s << "output_"
      << meshName;
    _context = s.str();

    PYLITH_METHOD_END;
} // open


// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriter::close(void) {
    _context = "";
    _isOpen = false;
} // close


// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriter::openTimeStep(const PylithScalar t,
                                         const topology::Mesh& mesh) {
    // Default: no implementation.
} // openTimeStep


// ----------------------------------------------------------------------
// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriter::closeTimeStep(void) {
    // Default: no implementation.
} // closeTimeStep


// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriter::DataWriter(const DataWriter& w) :
    _context(w._context),
    _isInfo(w._isInfo),
    _isOpen(w._isOpen) {}


// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::DataWriter::writePointNames(const pylith::string_vector& names,
                                            const topology::Mesh& mesh) {
    // Default: no implementation.

} // writePointNames


// ----------------------------------------------------------------------
// Create and populate PETSc global vector with coordinates of mesh vertices.
void
pylith::meshio::DataWriter::getCoordsGlobalVec(PetscVec* coordsGlobalVec,
                                               const pylith::topology::Mesh& mesh) {
    PYLITH_METHOD_BEGIN;

    assert(coordsGlobalVec);

    PetscDM dmMesh = mesh.getDM();assert(dmMesh);
    PetscDM dmCoord = NULL;

    PetscSection section = NULL;
    PetscSection newSection = NULL;
    PetscSection gsection = NULL;

    PetscSection subSection = NULL;
    PetscDMLabel subpointMap = NULL;
    PetscDMLabel subpointMapF = NULL;

    PylithInt dim, dimF, pStart, pEnd, qStart, qEnd, cStart, cEnd, cMax, vEnd, vMax = -1;

    PylithCallPetsc(DMGetCoordinateDM(dmMesh, &dmCoord));assert(dmCoord);
    PylithCallPetsc(DMPlexGetHeightStratum(dmCoord, 0, &cStart, &cEnd));
    PylithCallPetsc(DMPlexGetDepthStratum(dmCoord, 0, NULL, &vEnd));
    cMax = cStart;
    for (PetscInt cell = cStart; cell < cEnd; ++cell, ++cMax) {
        if (pylith::topology::MeshOps::isCohesiveCell(dmMesh, cell)) { break; }
    }
    PylithInt excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
    PylithInt numExcludes = (cMax < cEnd ? 1 : 0) + (vMax >= 0 ? 1 : 0);

    PylithCallPetsc(DMGetLocalSection(dmCoord, &section));
    PylithCallPetsc(DMGetDimension(dmMesh,  &dim));
    PylithCallPetsc(DMGetDimension(dmCoord, &dimF));
    PylithCallPetsc(DMPlexGetChart(dmMesh,  &pStart, &pEnd));
    PylithCallPetsc(DMPlexGetChart(dmCoord, &qStart, &qEnd));
    PylithCallPetsc(DMPlexGetSubpointMap(dmMesh,  &subpointMap));
    PylithCallPetsc(DMPlexGetSubpointMap(dmCoord, &subpointMapF));
    if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
        const PylithInt *indices = NULL;
        PetscIS subpointIS = NULL;
        PylithInt n = 0;

        PylithCallPetsc(PetscSectionGetChart(section, &qStart, &qEnd));
        PylithCallPetsc(DMPlexGetSubpointIS(dmMesh, &subpointIS));
        if (subpointIS) {
            PylithCallPetsc(ISGetLocalSize(subpointIS, &n));
            PylithCallPetsc(ISGetIndices(subpointIS, &indices));
        } // if
        PylithCallPetsc(PetscSectionCreate(mesh.getComm(), &subSection));
        PylithCallPetsc(PetscSectionSetChart(subSection, pStart, pEnd));
        for (PylithInt q = qStart; q < qEnd; ++q) {
            PylithInt dof, off, p;

            PylithCallPetsc(PetscSectionGetDof(section, q, &dof));
            if (dof) {
                PylithCallPetsc(PetscFindInt(q, n, indices, &p));
                if ((p >= pStart) && (p < pEnd)) {
                    PylithCallPetsc(PetscSectionSetDof(subSection, p, dof));
                    PylithCallPetsc(PetscSectionGetOffset(section, q, &off));
                    PylithCallPetsc(PetscSectionSetOffset(subSection, p, off));
                } // if
            } // if
        } // for
        if (subpointIS) {
            PylithCallPetsc(ISRestoreIndices(subpointIS, &indices));
        } // if
          /* No need to setup section */
        section = subSection;
        /* There are no excludes for surface meshes */
        numExcludes = 0;
    } // if

    PetscSF sf = NULL;
    PetscDM dmCoordGlobal = NULL;
    PylithCallPetsc(DMClone(dmCoord, &dmCoordGlobal));
    PylithCallPetsc(DMCopyDisc(dmCoord, dmCoordGlobal));
    PylithCallPetsc(PetscSectionClone(section, &newSection));
    PylithCallPetsc(DMSetLocalSection(dmCoordGlobal, newSection));
    PylithCallPetsc(PetscSectionDestroy(&newSection));
    PylithCallPetsc(DMGetPointSF(dmCoordGlobal, &sf));
    PylithCallPetsc(PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection));
    PylithCallPetsc(DMSetGlobalSection(dmCoordGlobal, gsection));
    PylithCallPetsc(PetscSectionDestroy(&gsection));

    PylithCallPetsc(VecDestroy(coordsGlobalVec));
    PylithCallPetsc(DMCreateGlobalVector(dmCoordGlobal, coordsGlobalVec));
    PylithCallPetsc(PetscObjectSetName((PetscObject) *coordsGlobalVec, "vertices"));

    PylithCallPetsc(PetscSectionDestroy(&subSection));

    InsertMode mode = INSERT_VALUES;
    PetscVec coordsLocalVec = NULL;
    PylithCallPetsc(DMGetCoordinatesLocal(dmMesh, &coordsLocalVec));
    PylithCallPetsc(DMLocalToGlobalBegin(dmCoordGlobal, coordsLocalVec, mode, *coordsGlobalVec));
    PylithCallPetsc(DMLocalToGlobalEnd(dmCoordGlobal, coordsLocalVec, mode, *coordsGlobalVec));

    PylithReal lengthScale = 1.0;
    PylithCallPetsc(DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale));
    PylithCallPetsc(VecScale(*coordsGlobalVec, lengthScale));
    PylithCallPetsc(DMDestroy(&dmCoordGlobal));

    PYLITH_METHOD_END;
} // getCoordsGlobalVec


// ----------------------------------------------------------------------
void
pylith::meshio::DataWriter::_writeVec(PetscVec vector,
                                      PetscViewer viewer) {
    // From plexhdf5.c DMPlexGlobalVectorView_HDF5_Internal

    /* To save vec in where we want, we create a new Vec (temp) with           */
    /* VecCreate(), wrap the vec data in temp, and call VecView(temp, viewer). */
    PetscVec temp;
    const PetscScalar *array;
    const char* vecName;
    PetscLayout map;

    PylithCallPetsc(VecCreate(PetscObjectComm((PetscObject)vector), &temp));
    PylithCallPetsc(PetscObjectGetName((PetscObject)vector, &vecName));
    PylithCallPetsc(PetscObjectSetName((PetscObject)temp, vecName));
    PylithCallPetsc(VecGetLayout(vector, &map));
    PylithCallPetsc(VecSetLayout(temp, map));
    PylithCallPetsc(VecSetUp(temp));
    PylithCallPetsc(VecGetArrayRead(vector, &array));
    PylithCallPetsc(VecPlaceArray(temp, array));
    PylithCallPetsc(VecView(temp, viewer));
    PylithCallPetsc(VecResetArray(temp));
    PylithCallPetsc(VecRestoreArrayRead(vector, &array));
    PylithCallPetsc(VecDestroy(&temp));
} // _writeVec


// End of file

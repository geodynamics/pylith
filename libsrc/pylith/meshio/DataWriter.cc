// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

#include <portinfo>

#include "DataWriter.hh" // Implementation of class methods

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
pylith::meshio::DataWriter::deallocate(void) {}


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
    PetscErrorCode err;

    err = DMGetCoordinateDM(dmMesh, &dmCoord);PYLITH_CHECK_ERROR(err);assert(dmCoord);
    err = DMPlexGetHeightStratum(dmCoord, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(dmCoord, 0, NULL, &vEnd);PYLITH_CHECK_ERROR(err);
    cMax = cStart;
    for (PetscInt cell = cStart; cell < cEnd; ++cell, ++cMax) {
        if (pylith::topology::MeshOps::isCohesiveCell(dmMesh, cell)) { break; }
    }
    PylithInt excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
    PylithInt numExcludes = (cMax < cEnd ? 1 : 0) + (vMax >= 0 ? 1 : 0);

    err = DMGetSection(dmCoord, &section);PYLITH_CHECK_ERROR(err);
    err = DMGetDimension(dmMesh,  &dim);PYLITH_CHECK_ERROR(err);
    err = DMGetDimension(dmCoord, &dimF);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(dmMesh,  &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(dmCoord, &qStart, &qEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(dmMesh,  &subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(dmCoord, &subpointMapF);PYLITH_CHECK_ERROR(err);
    if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
        const PylithInt *indices = NULL;
        PetscIS subpointIS = NULL;
        PylithInt n = 0;

        err = PetscSectionGetChart(section, &qStart, &qEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetSubpointIS(dmMesh, &subpointIS);PYLITH_CHECK_ERROR(err);
        if (subpointIS) {
            err = ISGetLocalSize(subpointIS, &n);PYLITH_CHECK_ERROR(err);
            err = ISGetIndices(subpointIS, &indices);PYLITH_CHECK_ERROR(err);
        } // if
        err = PetscSectionCreate(mesh.getComm(), &subSection);PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(subSection, pStart, pEnd);PYLITH_CHECK_ERROR(err);
        for (PylithInt q = qStart; q < qEnd; ++q) {
            PylithInt dof, off, p;

            err = PetscSectionGetDof(section, q, &dof);PYLITH_CHECK_ERROR(err);
            if (dof) {
                err = PetscFindInt(q, n, indices, &p);PYLITH_CHECK_ERROR(err);
                if ((p >= pStart) && (p < pEnd)) {
                    err = PetscSectionSetDof(subSection, p, dof);PYLITH_CHECK_ERROR(err);
                    err = PetscSectionGetOffset(section, q, &off);PYLITH_CHECK_ERROR(err);
                    err = PetscSectionSetOffset(subSection, p, off);PYLITH_CHECK_ERROR(err);
                } // if
            } // if
        } // for
        if (subpointIS) {
            err = ISRestoreIndices(subpointIS, &indices);PYLITH_CHECK_ERROR(err);
        } // if
          /* No need to setup section */
        section = subSection;
        /* There are no excludes for surface meshes */
        numExcludes = 0;
    } // if

    PetscSF sf = NULL;
    PetscDM dmCoordGlobal = NULL;
    err = DMClone(dmCoord, &dmCoordGlobal);PYLITH_CHECK_ERROR(err);
    err = DMCopyDisc(dmCoord, dmCoordGlobal);PYLITH_CHECK_ERROR(err);
    err = PetscSectionClone(section, &newSection);PYLITH_CHECK_ERROR(err);
    err = DMSetSection(dmCoordGlobal, newSection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&newSection);PYLITH_CHECK_ERROR(err);
    err = DMGetPointSF(dmCoordGlobal, &sf);PYLITH_CHECK_ERROR(err);
    err = PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection);PYLITH_CHECK_ERROR(err);
    err = DMSetGlobalSection(dmCoordGlobal, gsection);PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&gsection);PYLITH_CHECK_ERROR(err);

    err = VecDestroy(coordsGlobalVec);PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(dmCoordGlobal, coordsGlobalVec);PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) *coordsGlobalVec, "vertices");PYLITH_CHECK_ERROR(err);

    err = PetscSectionDestroy(&subSection);PYLITH_CHECK_ERROR(err);

    InsertMode mode = INSERT_VALUES;
    PetscVec coordsLocalVec = NULL;
    err = DMGetCoordinatesLocal(dmMesh, &coordsLocalVec);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalBegin(dmCoordGlobal, coordsLocalVec, mode, *coordsGlobalVec);PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(dmCoordGlobal, coordsLocalVec, mode, *coordsGlobalVec);PYLITH_CHECK_ERROR(err);

    PylithReal lengthScale = 1.0;
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = VecScale(*coordsGlobalVec, lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMDestroy(&dmCoordGlobal);PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // getCoordsGlobalVec


// End of file

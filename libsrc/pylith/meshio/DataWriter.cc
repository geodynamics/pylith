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

#include "DataWriter.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_BEGIN/END

// ----------------------------------------------------------------------
// Constructor
pylith::meshio::DataWriter::DataWriter(void) :
    _timeScale(1.0),
    _isInfo(false),
    _context("")
{ // constructor
} // constructor

// ----------------------------------------------------------------------
// Destructor
pylith::meshio::DataWriter::~DataWriter(void)
{ // destructor
    deallocate();
} // destructor

// ----------------------------------------------------------------------
// Deallocate PETSc and local data structures.
void
pylith::meshio::DataWriter::deallocate(void)
{ // deallocate
} // deallocate


// ----------------------------------------------------------------------
// Set time scale for simulation time.
void
pylith::meshio::DataWriter::timeScale(const PylithScalar value)
{ // timeScale
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
// Prepare for writing files.
void
pylith::meshio::DataWriter::open(const topology::Mesh& mesh,
                                 const bool isInfo,
                                 const char* label,
                                 const int labelId)
{ // open
    PYLITH_METHOD_BEGIN;

    _isInfo = isInfo;

    PetscDM dmMesh = mesh.dmMesh(); assert(dmMesh);
    const char* meshName = NULL;
    PetscObjectGetName((PetscObject) dmMesh, &meshName);

    std::ostringstream s;
    s << "output_"
      << meshName;
    if (label)
        s << "_" << label << labelId;
    _context = s.str();

    PYLITH_METHOD_END;
} // open

// ----------------------------------------------------------------------
// Close output files.
void
pylith::meshio::DataWriter::close(void)
{ // close
    _context = "";
} // close

// ----------------------------------------------------------------------
// Prepare file for data at a new time step.
void
pylith::meshio::DataWriter::openTimeStep(const PylithScalar t,
                                         const topology::Mesh& mesh,
                                         const char* label,
                                         const int labelId)
{ // openTimeStep
  // Default: no implementation.
} // openTimeStep

// ----------------------------------------------------------------------
// Cleanup after writing data for a time step.
void
pylith::meshio::DataWriter::closeTimeStep(void)
{ // closeTimeStep
  // Default: no implementation.
} // closeTimeStep

// ----------------------------------------------------------------------
// Copy constructor.
pylith::meshio::DataWriter::DataWriter(const DataWriter& w) :
    _isInfo(w._isInfo),
    _context(w._context)
{ // copy constructor
} // copy constructor


// ----------------------------------------------------------------------
// Write dataset with names of points to file.
void
pylith::meshio::DataWriter::writePointNames(const pylith::string_vector& names,
                                            const topology::Mesh& mesh)
{ // writePointNames
  // Default: no implementation.
} // writePointNames


// ----------------------------------------------------------------------
// Create and populate PETSc global vector with coordinates of mesh vertices.
void
pylith::meshio::DataWriter::getCoordsGlobalVec(PetscVec* coordsGlobalVec,
					       const pylith::topology::Mesh& mesh)
{ // getCoordsGlobalVec
    PYLITH_METHOD_BEGIN;
    
    assert(coordsGlobalVec);

    PetscDM dmMesh = mesh.dmMesh(); assert(dmMesh);
    PetscDM dmCoord = NULL;
    
    PetscSection section = NULL;
    PetscSection newSection = NULL;
    PetscSection gsection = NULL;
    
    PetscSection subSection = NULL;
    PetscDMLabel subpointMap = NULL;
    PetscDMLabel subpointMapF = NULL;
    
    PylithInt dim, dimF, pStart, pEnd, qStart, qEnd, cEnd, cMax, vEnd, vMax;
    PetscErrorCode err;
    
    err = DMGetCoordinateDM(dmMesh, &dmCoord); PYLITH_CHECK_ERROR(err); assert(dmCoord);
    err = DMPlexGetHeightStratum(dmCoord, 0, NULL, &cEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetDepthStratum(dmCoord, 0, NULL, &vEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(dmCoord, &cMax, NULL, NULL, &vMax); PYLITH_CHECK_ERROR(err);
    PylithInt excludeRanges[4] = {cMax, cEnd, vMax, vEnd};
    PylithInt numExcludes = (cMax >= 0 ? 1 : 0) + (vMax >= 0 ? 1 : 0);

    err = DMGetDefaultSection(dmCoord, &section); PYLITH_CHECK_ERROR(err);
    err = DMGetDimension(dmMesh,  &dim); PYLITH_CHECK_ERROR(err);
    err = DMGetDimension(dmCoord, &dimF); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(dmMesh,  &pStart, &pEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetChart(dmCoord, &qStart, &qEnd); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(dmMesh,  &subpointMap); PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSubpointMap(dmCoord, &subpointMapF); PYLITH_CHECK_ERROR(err);
    if (((dim != dimF) || ((pEnd-pStart) < (qEnd-qStart))) && subpointMap && !subpointMapF) {
        const PylithInt *indices = NULL;
        PetscIS subpointIS = NULL;
        PylithInt n = 0;

        err = PetscSectionGetChart(section, &qStart, &qEnd); PYLITH_CHECK_ERROR(err);
        err = DMPlexCreateSubpointIS(dmMesh, &subpointIS); PYLITH_CHECK_ERROR(err);
        if (subpointIS) {
            err = ISGetLocalSize(subpointIS, &n); PYLITH_CHECK_ERROR(err);
            err = ISGetIndices(subpointIS, &indices); PYLITH_CHECK_ERROR(err);
        } // if
        err = PetscSectionCreate(mesh.comm(), &subSection); PYLITH_CHECK_ERROR(err);
        err = PetscSectionSetChart(subSection, pStart, pEnd); PYLITH_CHECK_ERROR(err);
        for (PylithInt q = qStart; q < qEnd; ++q) {
            PylithInt dof, off, p;

            err = PetscSectionGetDof(section, q, &dof); PYLITH_CHECK_ERROR(err);
            if (dof) {
                err = PetscFindInt(q, n, indices, &p); PYLITH_CHECK_ERROR(err);
                if ((p >= pStart) && (p < pEnd)) {
                    err = PetscSectionSetDof(subSection, p, dof); PYLITH_CHECK_ERROR(err);
                    err = PetscSectionGetOffset(section, q, &off); PYLITH_CHECK_ERROR(err);
                    err = PetscSectionSetOffset(subSection, p, off); PYLITH_CHECK_ERROR(err);
                } // if
            } // if
        } // for
        if (subpointIS) {
            err = ISRestoreIndices(subpointIS, &indices); PYLITH_CHECK_ERROR(err);
            err = ISDestroy(&subpointIS); PYLITH_CHECK_ERROR(err);
        } // if
          /* No need to setup section */
        section = subSection;
        /* There are no excludes for surface meshes */
        numExcludes = 0;
    } // if

    PetscDS prob = NULL;
    PetscSF sf = NULL;
    PetscDM dmCoordGlobal = NULL;
    err = DMClone(dmCoord, &dmCoordGlobal); PYLITH_CHECK_ERROR(err);
    err = DMGetDS(dmCoord, &prob); PYLITH_CHECK_ERROR(err);
    err = DMSetDS(dmCoordGlobal, prob); PYLITH_CHECK_ERROR(err);
    err = PetscSectionClone(section, &newSection); PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultSection(dmCoordGlobal, newSection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&newSection); PYLITH_CHECK_ERROR(err);
    err = DMGetPointSF(dmCoordGlobal, &sf); PYLITH_CHECK_ERROR(err);
    err = PetscSectionCreateGlobalSectionCensored(section, sf, PETSC_TRUE, numExcludes, excludeRanges, &gsection);PYLITH_CHECK_ERROR(err);
    err = DMSetDefaultGlobalSection(dmCoordGlobal, gsection); PYLITH_CHECK_ERROR(err);
    err = PetscSectionDestroy(&gsection); PYLITH_CHECK_ERROR(err);

    err = VecDestroy(coordsGlobalVec); PYLITH_CHECK_ERROR(err);
    err = DMCreateGlobalVector(dmCoordGlobal, coordsGlobalVec); PYLITH_CHECK_ERROR(err);
    err = PetscObjectSetName((PetscObject) *coordsGlobalVec, "vertices"); PYLITH_CHECK_ERROR(err);

    err = PetscSectionDestroy(&subSection); PYLITH_CHECK_ERROR(err);

    
    InsertMode mode = INSERT_VALUES;
    PetscVec coordsLocalVec = NULL;
    err = DMGetCoordinatesLocal(dmMesh, &coordsLocalVec); PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalBegin(dmCoordGlobal, coordsLocalVec, mode, *coordsGlobalVec); PYLITH_CHECK_ERROR(err);
    err = DMLocalToGlobalEnd(dmCoordGlobal, coordsLocalVec, mode, *coordsGlobalVec); PYLITH_CHECK_ERROR(err);

    PylithReal lengthScale = 1.0;
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale); PYLITH_CHECK_ERROR(err);
    err = VecScale(*coordsGlobalVec, lengthScale); PYLITH_CHECK_ERROR(err);

    PYLITH_METHOD_END;
} // getCoordsGlobalVec

// End of file

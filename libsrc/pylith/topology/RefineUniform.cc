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

#include "pylith/topology/RefineUniform.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <cassert> // USES assert()
#include <sstream> // USES std::ostringstream
#include <stdexcept> // USES std::runtime_error

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class _RefineUniform {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt refine;
                static PylithInt refineFixCellLabel;
                static PylithInt refineCheckTopology;
            };

        }; // _RefineUniform
    } // topology
} // pylith

pylith::utils::EventLogger pylith::topology::_RefineUniform::Events::logger;
PylithInt pylith::topology::_RefineUniform::Events::refine;
PylithInt pylith::topology::_RefineUniform::Events::refineFixCellLabel;
PylithInt pylith::topology::_RefineUniform::Events::refineCheckTopology;

// ------------------------------------------------------------------------------------------------
void
pylith::topology::_RefineUniform::Events::init(void) {
    logger.setClassName("RefineUniform");
    logger.initialize();
    refine = logger.registerEvent("PL:RefineUniform:refine");
    refineFixCellLabel = logger.registerEvent("PL:RefineUniform:refineFixCellLabel");
    refineCheckTopology = logger.registerEvent("PL:RefineUniform:refineCheckTopology");
}


// ----------------------------------------------------------------------
// Constructor
pylith::topology::RefineUniform::RefineUniform(void) {
    _RefineUniform::Events::init();
}


// ----------------------------------------------------------------------
// Destructor
pylith::topology::RefineUniform::~RefineUniform(void) {
    deallocate();
}


// ----------------------------------------------------------------------
// Deallocate data structures.
void
pylith::topology::RefineUniform::deallocate(void) {
}


// ----------------------------------------------------------------------
// Refine mesh.
void
pylith::topology::RefineUniform::refine(Mesh* const newMesh,
                                        const Mesh& mesh,
                                        const int levels) {
    PYLITH_METHOD_BEGIN;
    _RefineUniform::Events::logger.eventBegin(_RefineUniform::Events::refine);

    if (levels < 1) {
        PYLITH_METHOD_END;
    } // if

    assert(newMesh);

    PetscErrorCode err;
    PetscDM dmOrig = mesh.getDM();assert(dmOrig);

    PetscInt meshDepth = 0;
    err = DMPlexGetDepth(dmOrig, &meshDepth);

    const int meshDim = mesh.getDimension();
    if (( meshDim > 0) && ( meshDepth != meshDim) ) {
        std::ostringstream msg;
        msg << "Mesh refinement for uninterpolated meshes not supported.\n"
            << "Turn on interpolated meshes using 'interpolate' mesh generator property.";
        throw std::runtime_error(msg.str());
    } // if

    // Refine, keeping original mesh intact.
    PetscDM dmNew = NULL;
    err = DMPlexSetRefinementUniform(dmOrig, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
    err = DMRefine(dmOrig, mesh.getComm(), &dmNew);PYLITH_CHECK_ERROR(err);

    for (int i = 1; i < levels; ++i) {
        PetscDM dmCur = dmNew;dmNew = NULL;
        err = DMPlexSetRefinementUniform(dmCur, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        err = DMRefine(dmCur, mesh.getComm(), &dmNew);PYLITH_CHECK_ERROR(err);

        err = DMDestroy(&dmCur);PYLITH_CHECK_ERROR(err);
    } // for
    err = DMPlexReorderCohesiveSupports(dmNew);PYLITH_CHECK_ERROR(err);

    newMesh->setDM(dmNew, "domain");

    _RefineUniform::Events::logger.eventBegin(_RefineUniform::Events::refineFixCellLabel);
    // Remove all non-cells from cells label
    const char* const labelName = pylith::topology::Mesh::cells_label_name;
    PetscDMLabel matidLabel = NULL;
    PetscIS valuesIS = NULL;
    const PetscInt *values = NULL;
    PetscInt cStart, cEnd, labelNumValues;
    err = DMPlexGetHeightStratum(dmNew, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMGetLabel(dmNew, labelName, &matidLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetNumValues(matidLabel, &labelNumValues);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetValueIS(matidLabel, &valuesIS);PYLITH_CHECK_ERROR(err);
    err = ISGetIndices(valuesIS, &values);PYLITH_CHECK_ERROR(err);
    for (PetscInt iValue = 0; iValue < labelNumValues; ++iValue) {
        PetscIS stratumIS = NULL;
        const PetscInt *points = NULL;
        const PetscInt value = values[iValue];
        PetscInt numPoints;
        err = DMLabelGetStratumSize(matidLabel, value, &numPoints);PYLITH_CHECK_ERROR(err);
        err = DMLabelGetStratumIS(matidLabel, value, &stratumIS);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(stratumIS, &points);PYLITH_CHECK_ERROR(err);
        for (PetscInt p = 0; p < numPoints; ++p) {
            const PetscInt point = points[p];
            if (( point < cStart) || ( point >= cEnd) ) {
                err = DMLabelClearValue(matidLabel, point, value);PYLITH_CHECK_ERROR(err);
            } // if
        } // for
        err = ISRestoreIndices(stratumIS, &points);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&stratumIS);PYLITH_CHECK_ERROR(err);
    } // for
    err = ISRestoreIndices(valuesIS, &values);PYLITH_CHECK_ERROR(err);
    err = ISDestroy(&valuesIS);PYLITH_CHECK_ERROR(err);
    _RefineUniform::Events::logger.eventEnd(_RefineUniform::Events::refineFixCellLabel);

    _RefineUniform::Events::logger.eventBegin(_RefineUniform::Events::refineCheckTopology);
    // Check consistency
    pylith::topology::MeshOps::checkTopology(*newMesh);
    _RefineUniform::Events::logger.eventEnd(_RefineUniform::Events::refineCheckTopology);

    // newMesh->view("REFINED_MESH", "::ascii_info_detail");

    _RefineUniform::Events::logger.eventEnd(_RefineUniform::Events::refine);
    PYLITH_METHOD_END;
} // refine


// End of file

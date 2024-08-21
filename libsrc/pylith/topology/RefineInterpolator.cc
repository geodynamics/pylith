// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

#include <portinfo>

#include "pylith/topology/RefineInterpolator.hh" // implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES MeshOps
#include "pylith/topology/Field.hh" // USES Field

#include "pylith/utils/error.hh" // USES PYLITH_METHOD_*
#include "pylith/utils/EventLogger.hh" // USES EventLogger

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
namespace pylith {
    namespace topology {
        class _RefineInterpolator {
public:

            class Events {
public:

                static
                void init(void);

                static pylith::utils::EventLogger logger;
                static PylithInt initialize;
                static PylithInt interpolate;
            };

        }; // _RefineInterpolator
    } // topology
} // pylith

pylith::utils::EventLogger pylith::topology::_RefineInterpolator::Events::logger;
PylithInt pylith::topology::_RefineInterpolator::Events::initialize;
PylithInt pylith::topology::_RefineInterpolator::Events::interpolate;

// ------------------------------------------------------------------------------------------------
void
pylith::topology::_RefineInterpolator::Events::init(void) {
    logger.setClassName("RefineInterpolator");
    logger.initialize();
    initialize = logger.registerEvent("PL:RefineInterpolator:initialize");
    interpolate = logger.registerEvent("PL:RefineInterpolator:interpolate");
}


// ------------------------------------------------------------------------------------------------
// Constructor
pylith::topology::RefineInterpolator::RefineInterpolator(void) {
    _RefineInterpolator::Events::init();
}


// ------------------------------------------------------------------------------------------------
// Destructor
pylith::topology::RefineInterpolator::~RefineInterpolator(void) {
    deallocate();
}


// ------------------------------------------------------------------------------------------------
// Deallocate data structures.
void
pylith::topology::RefineInterpolator::deallocate(void) {
    for (auto level : _levels) {
        DMDestroy(&level.dm);
        MatDestroy(&level.interpolateMatrix);
        VecDestroy(&level.vector);
    } // for
    _levels.resize(0);
} // deallocate


// ------------------------------------------------------------------------------------------------
// Get PETSc DM for output (finest level)
PetscDM
pylith::topology::RefineInterpolator::getOutputDM(void) {
    return (_levels.size() > 0) ? _levels[_levels.size()-1].dm : PETSC_NULLPTR;
}


// ------------------------------------------------------------------------------------------------
// Initialize interpolation to refined mesh.
void
pylith::topology::RefineInterpolator::initialize(const PetscDM& dmMesh,
                                                 const int refineLevels) {
    PYLITH_METHOD_BEGIN;
    _RefineInterpolator::Events::logger.eventBegin(_RefineInterpolator::Events::initialize);

    _levels.resize(refineLevels);
    PetscErrorCode err = PETSC_SUCCESS;

    PetscDM dmStart = pylith::topology::MeshOps::removeHangingCells(dmMesh);
    err = DMCopyDisc(dmMesh, dmStart);PYLITH_CHECK_ERROR(err);

    PetscDM dmPrev = dmStart;
    PetscReal lengthScale = 1.0;
    for (size_t iLevel = 0; iLevel < _levels.size(); ++iLevel) {
        _levels[iLevel].dm = PETSC_NULLPTR;
        _levels[iLevel].interpolateMatrix = PETSC_NULLPTR;
        _levels[iLevel].vector = PETSC_NULLPTR;

        err = DMPlexSetRefinementUniform(dmPrev, PETSC_TRUE);PYLITH_CHECK_ERROR(err);
        err = DMRefine(dmPrev, PetscObjectComm((PetscObject) dmMesh), &_levels[iLevel].dm);PYLITH_CHECK_ERROR(err);
        err = DMSetCoarseDM(_levels[iLevel].dm, dmPrev);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetScale(dmPrev, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
        err = DMPlexSetScale(_levels[iLevel].dm, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

#if 0 // needed for higher order coordinates (not needed for affine coordinates)
        PetscCall(DMPlexCreateCoordinateSpace(rdm, rd, PETSC_FALSE, NULL));
        PetscCall(PetscObjectSetName((PetscObject)rdm, "Refined Mesh with Linear Coordinates"));
        PetscCall(DMGetCoordinateDM(odm, &cdm));
        PetscCall(DMGetCoordinateDM(rdm, &rcdm));
        PetscCall(DMGetCoordinatesLocal(odm, &cl));
        PetscCall(DMGetCoordinatesLocal(rdm, &rcl));

#endif
        err = DMCopyDisc(dmPrev, _levels[iLevel].dm);PYLITH_CHECK_ERROR(err);
        err = DMGetGlobalVector(_levels[iLevel].dm, &_levels[iLevel].vector);PYLITH_CHECK_ERROR(err);
        err = PetscObjectReference((PetscObject) _levels[iLevel].vector);PYLITH_CHECK_ERROR(err);

        err = DMCreateInterpolation(dmPrev, _levels[iLevel].dm, &_levels[iLevel].interpolateMatrix, NULL);PYLITH_CHECK_ERROR(err);

        dmPrev = _levels[iLevel].dm;
    } // for

    _RefineInterpolator::Events::logger.eventEnd(_RefineInterpolator::Events::initialize);
    PYLITH_METHOD_END;
} // initialize


// ------------------------------------------------------------------------------------------------
// Interpolate field to fine mesh level.
void
pylith::topology::RefineInterpolator::interpolate(const PetscVec* vectorOut,
                                                  const PetscVec& vectorIn) {
    PYLITH_METHOD_BEGIN;
    _RefineInterpolator::Events::logger.eventBegin(_RefineInterpolator::Events::interpolate);
    assert(vectorOut);

    // Should be global vectors

    PetscVec vectorPrev = vectorIn;
    PetscErrorCode err = PETSC_SUCCESS;
    for (auto level : _levels) {
        err = MatMult(level.interpolateMatrix, vectorPrev, level.vector);PYLITH_CHECK_ERROR(err);
        vectorPrev = level.vector;
    } // for
    err = VecCopy(vectorPrev, *vectorOut);

    _RefineInterpolator::Events::logger.eventEnd(_RefineInterpolator::Events::interpolate);
    PYLITH_METHOD_END;
} // interpolate


// End of file

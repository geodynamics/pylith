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
#include "pylith/topology/FieldOps.hh" // USES FieldOps

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
// Get PETSc DM for input (coarsest level)
PetscDM
pylith::topology::RefineInterpolator::getInputDM(void) {
    PetscDM dmStart = PETSC_NULLPTR;
    if (_levels.size() > 0) {
        PylithCallPetsc(DMGetCoarseDM(_levels[0].dm, &dmStart));
    } // if
    return dmStart;
}


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
                                                 const int refineLevels,
                                                 const int outputBasisOrder,
                                                 const pylith::topology::FieldBase::Description& description,
                                                 const pylith::topology::FieldBase::Discretization& discretization) {
    PYLITH_METHOD_BEGIN;
    _RefineInterpolator::Events::logger.eventBegin(_RefineInterpolator::Events::initialize);

    _levels.resize(refineLevels);

    PetscDM dmStart = pylith::topology::MeshOps::removeHangingCells(dmMesh);assert(dmStart);
    PylithCallPetsc(DMCopyDisc(dmMesh, dmStart));

    PetscDM dmPrev = dmStart;
    PetscReal lengthScale = 1.0;
    MPI_Comm comm = PetscObjectComm((PetscObject) dmMesh);
    for (size_t iLevel = 0; iLevel < _levels.size(); ++iLevel) {
        _levels[iLevel].dm = PETSC_NULLPTR;
        _levels[iLevel].interpolateMatrix = PETSC_NULLPTR;
        _levels[iLevel].vector = PETSC_NULLPTR;

        PylithCallPetsc(DMPlexSetRefinementUniform(dmPrev, PETSC_TRUE));
        PylithCallPetsc(DMRefine(dmPrev, comm, &_levels[iLevel].dm));
        PylithCallPetsc(DMSetCoarseDM(_levels[iLevel].dm, dmPrev));
        PylithCallPetsc(DMPlexGetScale(dmPrev, PETSC_UNIT_LENGTH, &lengthScale));
        PylithCallPetsc(DMPlexSetScale(_levels[iLevel].dm, PETSC_UNIT_LENGTH, lengthScale));
        PylithCallPetsc(DMPlexReorderSetDefault(_levels[iLevel].dm, DM_REORDER_DEFAULT_FALSE));

#if 0 // needed for higher order coordinates (not needed for affine coordinates)
        PetscCall(DMPlexCreateCoordinateSpace(rdm, rd, PETSC_FALSE, NULL));
        PetscCall(PetscObjectSetName((PetscObject)rdm, "Refined Mesh with Linear Coordinates"));
        PetscCall(DMGetCoordinateDM(odm, &cdm));
        PetscCall(DMGetCoordinateDM(rdm, &rcdm));
        PetscCall(DMGetCoordinatesLocal(odm, &cl));
        PetscCall(DMGetCoordinatesLocal(rdm, &rcl));
#endif

        if (iLevel < _levels.size()-1) {
            PylithCallPetsc(DMCopyDisc(dmPrev, _levels[iLevel].dm));
        } else {
            const PetscInt minBasisOrder = PETSC_DETERMINE;
            const PetscInt maxBasisOrder = outputBasisOrder;
            PylithCallPetsc(DMCopyFields(dmPrev, minBasisOrder, maxBasisOrder, _levels[iLevel].dm));
            PylithCallPetsc(DMCopyDS(dmPrev, minBasisOrder, maxBasisOrder, _levels[iLevel].dm));
        } // else
        PylithCallPetsc(DMCreateGlobalVector(_levels[iLevel].dm, &_levels[iLevel].vector));
        PylithCallPetsc(DMCreateInterpolation(dmPrev, _levels[iLevel].dm, &_levels[iLevel].interpolateMatrix, NULL));

        dmPrev = _levels[iLevel].dm;
    } // for
    PylithCallPetsc(DMDestroy(&dmStart));

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

    PetscVec vectorPrev = vectorIn;
    for (auto level : _levels) {
        PylithCallPetsc(MatMult(level.interpolateMatrix, vectorPrev, level.vector));
        vectorPrev = level.vector;
    } // for
    PylithCallPetsc(VecCopy(vectorPrev, *vectorOut));

    _RefineInterpolator::Events::logger.eventEnd(_RefineInterpolator::Events::interpolate);
    PYLITH_METHOD_END;
} // interpolate


// End of file

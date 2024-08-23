// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/utils/petscfwd.h" // USES PetscDM

#include "pylith/topology/FieldBase.hh" // USES FieldBase

#include <vector> // HASA std::vector

/// @brief  Interpolate fields to uniformly refined mesh.
class pylith::topology::RefineInterpolator {
    friend class TestRefineUniform; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    RefineInterpolator(void);

    /// Destructor
    ~RefineInterpolator(void);

    /// Deallocate data structures.
    void deallocate(void);

    /** Get PETSc DM for output (finest level)
     *
     */
    PetscDM getOutputDM(void);

    /** Initialize interpolation to refined mesh.
     *
     * @param[in] dmMesh PETSc DM for starting point of refinement.
     * @param[in] refineLevels Number of levels of mesh refinement.
     * @param[in] outputBasisOrder Basis order for output.
     */
    void initialize(const PetscDM& dmMesh,
                    const int refineLevels,
                    const int outputBasisOrder,
                    const pylith::topology::FieldBase::Description& description,
                    const pylith::topology::FieldBase::Discretization& discretization);

    /** Interpolate field to fine mesh level.
     *
     * @param[out] vectorOut PETSc Vec with interpolated field.
     * @param[in] vectorIn PETSc Vec with field to interpolate.
     */
    void interpolate(const PetscVec* vectorOut,
                     const PetscVec& vectorIn);

    // PRIvATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    struct Level {
        PetscDM dm;
        PetscMat interpolateMatrix;
        PetscVec vector;
    };

    std::vector<Level> _levels; ///< Information at each refinement level.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    RefineInterpolator(const RefineInterpolator&); ///< Not implemented
    const RefineInterpolator& operator=(const RefineInterpolator&); ///< Not implemented

}; // RefineInterpolator

// End of file

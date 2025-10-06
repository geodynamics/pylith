// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/utils/utilsfwd.hh" // forward declarations

#include "pylith/utils/petscfwd.h"

#include "petscts.h"

class pylith::utils::TSAdaptImpulse {
    friend class TestTSAdapImpulse; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /** Set time step for impulse portion.
     *
     * @param[in] timestep Time step (nondimensional).
     */
    static
    void setImpulseTimeStep(const double timestep);

    /** Get time step for impulse portion.
     *
     * @returns Time step (nondimensional).
     */
    static
    double getImpulseTimeStep(void);

    /** Set PETSc TS adapt.
     *
     * @param[inout] ts PETSc TS.
     */
    static
    void set(PetscTS ts);

    /** PETSc adaptive time stepper.
     *
     */
    static
    PetscErrorCode TSAdaptChoose(TSAdapt adapt,
                                 PetscTS ts,
                                 PetscReal h,
                                 PetscInt *next_sc,
                                 PetscReal *next_h,
                                 PetscBool *accept,
                                 PetscReal *wlte,
                                 PetscReal *wltea,
                                 PetscReal *wlter);

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    static double _impulseStep; ///< Time step for impulse portion.

    // Not impemented /////////////////////////////////////////////////////////////////////////////
private:

    TSAdaptImpulse(void); ///< Not implemented.
};

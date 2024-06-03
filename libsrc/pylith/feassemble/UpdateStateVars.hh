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

#include "pylith/feassemble/feassemblefwd.hh"

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/petscfwd.h" // USES PetscIS, PetscDM, PetscVec

class pylith::feassemble::UpdateStateVars : public pylith::utils::GenericComponent {
    friend class TestUpdateStateVars; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    UpdateStateVars(void);

    /// Destructor.
    virtual ~UpdateStateVars(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Get PETSc DM associated with state variables.
     *
     * @returns PETSc DM for state variables.
     */
    PetscDM stateVarsDM(void);

    /** Get PETSc local vector associated with state variables.
     *
     * @returns PETSc local vector with state variables.
     */
    PetscVec stateVarsLocalVector(void);

    /** Initialize layout for updating state variables.
     *
     * @param[in] auxiliaryField Auxiliary field containing state variables.
     */
    void initialize(const pylith::topology::Field& auxiliaryField);

    /** Extract current state variables in auxiliary field in preparation for computing new ones.
     *
     * @param[inout] auxiliaryField Auxiliary field containing state variables.
     */
    void prepare(pylith::topology::Field* auxiliaryField);

    /** Update state variables in auxiliary field after computing them.
     *
     * @param[inout] auxiliaryField Auxiliary field containing state variables.
     */
    void restore(pylith::topology::Field* auxiliaryField);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PetscIS _stateVarsIS; ///< Petsc IS for state vars in auxiliary field.
    PetscDM _stateVarsDM; ///< Petsc DM for state vars subfield.
    PetscVec _stateVarsVecLocal; ///< Petsc Vec with global vector for state vars.
    PetscVec _stateVarsVecGlobal; ///< Petsc Vec with global vector for state vars.
    PetscVec _auxiliaryFieldVecGlobal; ///< Petsc Vec with global vector for auxiliary field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    UpdateStateVars(const UpdateStateVars &); ///< Not implemented
    const UpdateStateVars& operator=(const UpdateStateVars&); ///< Not implemented

}; // class UpdateStateVars

// End of file

// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// ----------------------------------------------------------------------
//

/** @file libsrc/feassemble/UpdateStateVars.hh
 *
 * @brief C++ object for managing updating state variables.
 */

#if !defined(pylith_feassemble_updatestatevars_hh)
#define pylith_feassemble_updatestatevars_hh

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

#endif // pylith_feassemble_updatestatevars_hh

// End of file

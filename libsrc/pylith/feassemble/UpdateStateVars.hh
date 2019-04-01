// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
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

    /** Initialize layout for updating state variables.
     *
     * @param[in] auxiliaryField Auxiliary field containing state variables.
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& auxiliaryField,
                    const pylith::topology::Field& solution);

    /** Setup values for updating state variables.
     *
     * @param[inout] auxiliaryField Auxiliary field containing state variables.
     * @param[in] solution Solution field.
     */
    void prepareValues(pylith::topology::Field* auxiliaryField,
                       const pylith::topology::Field& solution);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PetscIS* superIS; ///< Petsc IS for ??
    PetscDM superDM; ///< Petsc DM for ??
    PetscIS stateVarIS; ///< Petsc IS for state vars in auxiliary field.
    PetscDM stateVarDM; ///< Petsc DM for state vars subfield.
    PetscVec stateVarsSolnVecLocal; ///< Petsc Vec with local vector for state vars and solution field.
    PetscVec stateVarsSolnVecGlobal; ///< Petsc Vec with global vector for state vars and solution field.
    PetscVec stateVarsVecGlobal; ///< Petsc Vec with global vector for state vars.
    PetscVec auxiliaryFieldVecGlobal; ///< Petsc Vec with global vector for auxiliary field.
    PetscVec solutionVecGlobal; ///< Petsc Vec with global vector for solution field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    UpdateStateVars(const UpdateStateVars &); ///< Not implemented
    const UpdateStateVars& operator=(const UpdateStateVars&); ///< Not implemented

}; // class UpdateStateVars

#endif // pylith_feassemble_updatestatevars_hh

// End of file

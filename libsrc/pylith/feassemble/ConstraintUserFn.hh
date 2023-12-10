// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::ConstraintUserFn : public pylith::feassemble::Constraint {
    friend class TestConstraintUserFn; // unit testing
    friend class _ConstraintUserFn; /// private utility class

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] physics Physics implemented by constraint.
     */
    ConstraintUserFn(pylith::problems::Physics* const physics);

    /// Destructor.
    virtual ~ConstraintUserFn(void);

    /** Set user function specifying constrained values.
     *
     * @param[in] fn Function specifying contrained values.
     */
    void setUserFn(const PetscUserFieldFunc fn);

    /** Set user function time derivative specifying constrained values.
     *
     * @param[in] fnDot Function specifying contrained values time derivative.
     */
    void setUserFnDot(const PetscUserFieldFunc fnDot);

    /** Initialize constraint.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution);

    /** Set constrained values in solution field.
     *
     * @param[inout] integrationData Data needed to integrate governing equation.
     */
    void setSolution(pylith::feassemble::IntegrationData* integrationData);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    PetscUserFieldFunc _fn; ///< Function for computing constrained values.
    PetscUserFieldFunc _fnDot; ///< Function for computing constrained values.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ConstraintUserFn(void); /// Not implemented.
    ConstraintUserFn(const Constraint &); ///< Not implemented
    const ConstraintUserFn& operator=(const ConstraintUserFn&); ///< Not implemented

}; // class ConstraintUserFn

// End of file

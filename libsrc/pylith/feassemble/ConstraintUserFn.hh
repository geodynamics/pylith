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

/** @file libsrc/feassemble/Constraint.hh
 *
 * @brief C++ class for constraining degrees of freedom in the solution via a user-specified analytical function.
 */

#if !defined(pylith_feassemble_constraintuserfn_hh)
#define pylith_feassemble_constraintuserfn_hh

#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::ConstraintUserFn : public pylith::feassemble::Constraint {
    friend class TestConstraintUserFn; // unit testing

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
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    virtual
    void setSolution(pylith::topology::Field* solution,
                     const double t);

    /** Set constrained values time derivative in solution field.
     *
     * @param[out] solutionDot Solution field.
     * @param[in] t Current time.
     */
    virtual
    void setSolutionDot(pylith::topology::Field* solutionDot,
                        const double t);

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

#endif // pylith_feassemble_constraintuserfn_hh

// End of file

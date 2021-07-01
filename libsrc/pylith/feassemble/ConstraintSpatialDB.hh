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
 * @brief C++ class for constraining degrees of freedom in the solution via an auxiliary field constructed from
 *  spatial database(s).
 */

#if !defined(pylith_feassemble_constraintspatialdb_hh)
#define pylith_feassemble_constraintspatialdb_hh

#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/problems/problemsfwd.hh" // HASA Physics
#include "pylith/topology/topologyfwd.hh" // USES Field

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

class pylith::feassemble::ConstraintSpatialDB : public pylith::feassemble::Constraint {
    friend class TestConstraintSpatialDB; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] physics Physics implemented by constraint.
     */
    ConstraintSpatialDB(pylith::problems::Physics* const physics);

    /// Destructor.
    ~ConstraintSpatialDB(void);

    /** Set constraint kernel.
     *
     * @param kernel Kernel to compute constrained value from auxiliary field.
     */
    void setKernelConstraint(const PetscBdPointFunc kernel);

    /** Initialize constraint.
     *
     * @param[in] solution Solution field (layout).
     */
    void initialize(const pylith::topology::Field& solution);

    /** Update auxiliary field values to current time.
     *
     * @param[in] t Current time.
     */
    void updateState(const PylithReal t);

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    void setSolution(pylith::topology::Field* solution,
                     const double t);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Set constants used in finite-element kernels.
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     */
    void _setKernelConstants(const pylith::topology::Field& solution,
                             const PylithReal dt) const;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PetscBdPointFunc _kernelConstraint; ///< Kernel for computing constrained values from auxiliary field.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ConstraintSpatialDB(void); /// Not implemented.
    ConstraintSpatialDB(const ConstraintSpatialDB &); ///< Not implemented
    const ConstraintSpatialDB& operator=(const ConstraintSpatialDB&); ///< Not implemented

}; // class Constraint

#endif // pylith_feassemble_constraintspatialdb_hh

// End of file

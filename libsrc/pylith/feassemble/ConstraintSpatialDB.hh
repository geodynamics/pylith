// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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

    /** Set auxiliary field values for current time.
     *
     * @param[in] t Current time.
     */
    void setState(const PylithReal t);

    /** Set constrained values in solution field.
     *
     * @param[inout] integrationData Data needed to integrate governing equation.
     */
    void setSolution(pylith::feassemble::IntegrationData* integrationData);

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

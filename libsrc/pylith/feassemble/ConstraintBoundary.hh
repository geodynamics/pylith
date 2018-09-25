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

/** @file libsrc/feassemble/ConstraintBoundary.hh
 *
 * @brief C++ implementation of constraining degrees of freedom on external boundary.
 */

#if !defined(pylith_feassemble_constraintboundary_hh)
#define pylith_feassemble_constraintboundary_hh

#include "feassemblefwd.hh" // forward declarations

#include "pylith/feassemble/Constraint.hh" // ISA Constraint

class pylith::feassemble::ConstraintBoundary : public pylith::feassemble::Constraint {
    friend class TestConstraintBoundary; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] physics Physics implemented by constraint.
     */
    ConstraintBoundary(pylith::problems::Physics* const physics);

    /// Destructor.
    ~ConstraintBoundary(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Get mesh associated with constrained boundary.
     *
     * @returns Mesh associated with constrained boundary.
     */
    const pylith::topology::Mesh& getConstraintDomainMesh(void) const;

    /** Initialize boundary condition.
     *
     * @param[in] solution Solution field.
     */
    void initialize(const pylith::topology::Field& solution);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _boundaryMesh; ///< Boundary mesh.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ConstraintBoundary(void); ///< Not implemented.
    ConstraintBoundary(const ConstraintBoundary&); ///< Not implemented.
    const ConstraintBoundary& operator=(const ConstraintBoundary&); ///< Not implemented.

};

// class ConstraintBoundary

#endif // pylith_feassemble_constraintboundary_hh

// End of file

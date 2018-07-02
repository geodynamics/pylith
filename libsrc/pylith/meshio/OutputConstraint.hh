// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputConstraint.hh
 *
 * @brief C++ object for managing output of the solution over the domain of a constraint.
 */

#if !defined(pylith_meshio_outputconstraint_hh)
#define pylith_meshio_outputconstraint_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "OutputManager.hh" // ISA OutputManager

#include "pylith/utils/array.hh" // HASA string_vector

// OutputConstraint -----------------------------------------------------
/** @brief C++ object for managing output of the solution over the domain of a constraint.
 */
class pylith::meshio::OutputConstraint : public pylith::meshio::OutputManager {
    friend class TestOutputConstraint;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] constraint Constraint to observe.
     */
    OutputConstraint(pylith::feassemble::ConstraintPointwise* const constraint);

    /// Destructor
    ~OutputConstraint(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Verify configuration.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    /// Write solution at time step.
    void _writeInfo(void);

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void _writeDataStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    pylith::feassemble::ConstraintPointwise* const _constraint;

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputConstraint(const OutputConstraint&);   ///< Not implemented.
    const OutputConstraint& operator=(const OutputConstraint&);   ///< Not implemented

}; // OutputConstraint

#endif // pylith_meshio_outputconstraint_hh

// End of file

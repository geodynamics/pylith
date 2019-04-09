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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/InitialCondition.hh
 *
 * @brief C++ abstract base class for specifying initial conditions.
 */
#if !defined(pylith_problems_initialcondition_hh)
#define pylith_problems_initialcondition_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/utils/arrayfwd.hh" // USES string_vector

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

class pylith::problems::InitialCondition : public pylith::utils::PyreComponent {
    friend class TestInitialCondition; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    InitialCondition(void);

    /// Destructor
    virtual ~InitialCondition(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set fields for initial condition.
     *
     * @param[in] fields Array of names of fields.
     * @param[in] numFields Number of fields.
     */
    void setFields(const char* fields[],
                   const int numFields);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Set solution to values for initial condition.
     *
     * @param[out] solution Solution field.
     * @param[in] normalizer Nondimensionalization.
     */
    virtual
    void setValues(pylith::topology::Field* solution,
                   const spatialdata::units::Nondimensional& normalizer) = 0;

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::string_vector _fields; ///< Names of fields for initial conditions.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    InitialCondition(const InitialCondition&); ///< Not implemented
    const InitialCondition& operator=(const InitialCondition&); ///< Not implemented

}; // InitialCondition

#endif // pylith_problems_initialcondition_hh

// End of file

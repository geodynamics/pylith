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
 * @file libsrc/problems/InitialConditions.hh
 *
 * @brief C++ abstract base class for specifying initial conditions.
 */
#if !defined(pylith_problems_initialconditions_hh)
#define pylith_problems_initialconditions_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

class pylith::problems::InitialConditions : public pylith::utils::PyreComponent {
    friend class TestInitialConditions; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    InitialConditions(void);

    /// Destructor
    virtual ~InitialConditions(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set solver type.
     *
     * @param[out] solution Solution field.
     * @param[in] normalizer Nondimensionalization.
     */
    virtual
    void setValues(pylith::topology::Field* solution,
                   const spatialdata::units::Nondimensional& normalizer) = 0;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    InitialConditions(const InitialConditions&); ///< Not implemented
    const InitialConditions& operator=(const InitialConditions&); ///< Not implemented

}; // InitialConditions

#endif // pylith_problems_initialconditions_hh

// End of file

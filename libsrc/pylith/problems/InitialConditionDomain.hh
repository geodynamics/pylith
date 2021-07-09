// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/problems/InitialConditionDomain.hh
 *
 * @brief C++ object for specifying initial conditions over the entire domain.
 */
#if !defined(pylith_problems_initialconditionsdomain_hh)
#define pylith_problems_initialconditionsdomain_hh

#include "InitialCondition.hh" // ISA InitialCondition

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB

class pylith::problems::InitialConditionDomain : public pylith::problems::InitialCondition {
    friend class TestInitialConditionDomain; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    InitialConditionDomain(void);

    /// Destructor
    virtual ~InitialConditionDomain(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set spatial database holding initial conditions.
     *
     * @param[in] db Spatial database holding initial conditions.
     */
    void setDB(spatialdata::spatialdb::SpatialDB* db);

    /** Set solution to values for initial condition.
     *
     * @param[out] solution Solution field.
     * @param[in] normalizer Nondimensionalization.
     */
    void setValues(pylith::topology::Field* solution,
                   const spatialdata::units::Nondimensional& normalizer);

    // PRIVATE MEMEBRS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database with values for initial condition.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    InitialConditionDomain(const InitialConditionDomain&); ///< Not implemented
    const InitialConditionDomain& operator=(const InitialConditionDomain&); ///< Not implemented

}; // InitialConditionDomain

#endif // pylith_problems_initialconditionsdomain_hh

// End of file

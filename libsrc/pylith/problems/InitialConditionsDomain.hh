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
 * @file libsrc/problems/InitialConditionsDomain.hh
 *
 * @brief C++ object for specifying initial conditions over the entire domain.
 */
#if !defined(pylith_problems_initialconditionsdomain_hh)
#define pylith_problems_initialconditionsdomain_hh

#include "InitialConditions.hh" // ISA InitialConditions

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB

class pylith::problems::InitialConditionsDomain : public pylith::problems::InitialConditions {
    friend class TestInitialConditionsDomain; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    InitialConditionsDomain(void);

    /// Destructor
    virtual ~InitialConditionsDomain(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set spatial database holding initial conditions.
     *
     * @param[in] db Spatial database holding initial conditions.
     */
    void setDB(spatialdata::spatialdb::SpatialDB* db);

    /** Set solver type.
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

    InitialConditionsDomain(const InitialConditionsDomain&); ///< Not implemented
    const InitialConditionsDomain& operator=(const InitialConditionsDomain&); ///< Not implemented

}; // InitialConditionsDomain

#endif // pylith_problems_initialconditionsdomain_hh

// End of file

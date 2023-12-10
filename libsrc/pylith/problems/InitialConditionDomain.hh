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

#include "pylith/problems/InitialCondition.hh" // ISA InitialCondition

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

// End of file

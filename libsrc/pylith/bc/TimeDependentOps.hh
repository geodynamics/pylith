// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2026, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/bc/bcfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES pylith::real

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES TimeHistoryDB

class pylith::bc::TimeDependentOps {
    friend class TestTimeDependentOps; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Update auxiliary field using time history database.
     *
     * @param[inout] auxiliaryField Auxiliary field to update.
     * @param[in] t Current nondimensional time.
     * @param[in] timeScale Scale for dimensionalizing time.
     * @param[in] dbTimeHistory Time history database.
     */
    static
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const pylith::real t,
                              const pylith::real timeScale,
                              const std::shared_ptr<spatialdata::spatialdb::TimeHistory>& dbTimeHistory);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    TimeDependentOps(void) = delete;
    TimeDependentOps(const TimeDependentOps&) = delete;
    const TimeDependentOps& operator=(const TimeDependentOps&) = delete;

}; // TimeDependentOps

// End of file

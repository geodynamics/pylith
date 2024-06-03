// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2024, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/faults/faultsfwd.hh" // forward declarations

#include "pylith/feassemble/feassemblefwd.hh" // USES IntegrationData
#include "pylith/topology/topologyfwd.hh" // USES Field

class pylith::faults::FaultOps {
    friend class TestFaultOps; // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Create weighting vector for dynamic prescribed slip integration.
     *
     * @param[in] integrationData Data for finite-element integration.
     */
    static
    void createDAEMassWeighting(pylith::feassemble::IntegrationData* integrationData);

    /** Update weighting vector for dynamic prescribed slip integration.
     *
     * @param[in] integrationData Data for finite-element integration.
     */
    static
    void updateDAEMassWeighting(pylith::feassemble::IntegrationData* integrationData);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    FaultOps(void); ///< Not implemented.
    FaultOps(const FaultOps&); ///< Not implemented.
    const FaultOps& operator=(const FaultOps&); ///< Not implemented.

}; // FaultOps

// End of file

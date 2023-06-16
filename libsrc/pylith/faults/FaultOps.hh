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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/faults/FaultOps.hh
 *
 * @brief Operations related to faults.
 */

#if !defined(pylith_faults_faultops_hh)
#define pylith_faults_faultops_hh

// Include directives ---------------------------------------------------
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

#endif // pylith_faults_fieldOps_hh

// End of file

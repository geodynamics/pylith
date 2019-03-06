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
 * @file libsrc/problems/ObserverSolnStub.hh
 *
 * @brief Minimal C++ implementation of ObserverSoln to allow unit tests using ObserverSoln objects.
 */

#if !defined(pylith_feassemble_observersolnstub_hh)
#define pylith_feassemble_observersolnstub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/problems/ObserverSoln.hh" // ISA ObserverSoln

class pylith::problems::ObserverSolnStub : public pylith::problems::ObserverSoln {
    friend class TestObserverSolnStub; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserverSolnStub(void);

    /// Destructor
    virtual ~ObserverSolnStub(void);

    /** Verify observer is compatible with solution.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

    /** Receive update (subject of observer).
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserverSolnStub(const ObserverSolnStub&); ///< Not implemented.
    const ObserverSolnStub& operator=(const ObserverSolnStub&); ///< Not implemented

}; // ObserverSolnStub

#endif // pylith_feassemble_observersolnstub_hh

// End of file

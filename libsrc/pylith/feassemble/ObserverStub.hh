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
 * @file libsrc/feassemble/ObserverStub.hh
 *
 * @brief Minimal C++ implementation of Observer to allow unit tests using Observer objects.
 */

#if !defined(pylith_feassemble_observerstub_hh)
#define pylith_feassemble_observerstub_hh

#include "pylith/feassemble/Observer.hh" // ISA Observer

class pylith::feassemble::ObserverStub : public pylith::feassemble::Observer {
    friend class TestObserverStub; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserverStub(void);

    /// Destructor
    virtual ~ObserverStub(void);

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
     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
     */
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution,
                const bool infoOnly=false);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserverStub(const ObserverStub&); ///< Not implemented.
    const ObserverStub& operator=(const ObserverStub&); ///< Not implemented

}; // ObserverStub

class pylith::feassemble::ObserverStubException {
    // PUBLIC ENUMS ////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum MethodEnum {
        VERIFIED=0,
        UPDATED=1,
    }; // MethodEnum

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Constructor.
     *
     * @param[in] value Method called.
     */
    ObserverStubException(const MethodEnum value);

    /// Destructor
    ~ObserverStubException(void);

    /** Get method called.
     *
     * @returns Method called.
     */
    MethodEnum getMethodCalled(void) const;

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
private:

    MethodEnum _methodCalled;
}; // ObserverStubException

#endif // pylith_feassemble_observerstub_hh

// End of file

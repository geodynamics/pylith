// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

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
    ~ObserverSolnStub(void);

    /** Set time scale.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    void setTimeScale(const PylithReal value);

    /** Get time scale.
     *
     * @returns Time scale for dimensionalizing time.
     */
    PylithReal getTimeScale(void) const;

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
     * @param[in] notification Type of notification.
     */
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution,
                const NotificationType notification);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _timeScale; ///< Time scale.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserverSolnStub(const ObserverSolnStub&); ///< Not implemented.
    const ObserverSolnStub& operator=(const ObserverSolnStub&); ///< Not implemented

}; // ObserverSolnStub

#endif // pylith_feassemble_observersolnstub_hh

// End of file

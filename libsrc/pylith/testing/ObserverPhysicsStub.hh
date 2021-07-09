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
 * @file libsrc/problems/ObserverPhysicsStub.hh
 *
 * @brief Minimal C++ implementation of Observer to allow unit tests using Observer objects.
 */

#if !defined(pylith_problems_observerphysicsstub_hh)
#define pylith_problems_observerphysicsstub_hh

#include "pylith/testing/testingfwd.hh" // forward declarations

#include "pylith/problems/ObserverPhysics.hh" // ISA ObserverPhysics

class pylith::problems::ObserverPhysicsStub : public pylith::problems::ObserverPhysics {
    friend class TestObserverPhysicsStub; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserverPhysicsStub(void);

    /// Destructor
    ~ObserverPhysicsStub(void);

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
     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
     */
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution,
                const bool infoOnly);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    PylithReal _timeScale; ///< Time scale.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserverPhysicsStub(const ObserverPhysicsStub&); ///< Not implemented.
    const ObserverPhysicsStub& operator=(const ObserverPhysicsStub&); ///< Not implemented

}; // ObserverPhysicsStub

#endif // pylith_problems_observerphysicsstub_hh

// End of file

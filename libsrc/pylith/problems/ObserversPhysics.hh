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
 * @file libsrc/problems/ObserversPhysics.hh
 *
 * @brief Collection of observers for an object.
 */

#if !defined(pylith_problems_observersphysics_hh)
#define pylith_problems_observersphysics_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/feassemble/feassemblefwd.hh" // USES PhysicsImplementation
#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES PylithReal, PylithInt

#include <set> // USES std::set

class pylith::problems::ObserversPhysics : public pylith::utils::GenericComponent {
    friend class TestObserversPhysics; // unit testing
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserversPhysics(void);

    /// Destructor
    virtual ~ObserversPhysics(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Register observer to receive notifications.
     *
     * ObserversPhysics are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(pylith::problems::ObserverPhysics* observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(pylith::problems::ObserverPhysics* observer);

    /** Get number of observers.
     *
     * @returns Number of observers.
     */
    size_t size(void) const;

    /** Set physics implementation in observers (for callbacks)
     *
     * @param[in] physics Physics implementation being observed.
     */
    void setPhysicsImplementation(const pylith::feassemble::PhysicsImplementation* const physics);

    /** Set time scale in observers.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    void setTimeScale(const PylithReal value);

    /** Verify observers are compatible.
     *
     * @param[in] solution Solution field.
     */
    void verifyObservers(const pylith::topology::Field& solution) const;

    /** Send observers an update.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
     */
    void notifyObservers(const PylithReal t,
                         const PylithInt tindex,
                         const pylith::topology::Field& solution,
                         const bool infoOnly);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    typedef std::set<pylith::problems::ObserverPhysics*>::iterator iterator; ///< Iterator.
    std::set<pylith::problems::ObserverPhysics*> _observers; ///< Subscribers of updates.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserversPhysics(const ObserversPhysics&); ///< Not implemented.
    const ObserversPhysics& operator=(const ObserversPhysics&); ///< Not implemented

};

// ObserversPhysics

#endif // pylith_problems_observersphysics_hh

// End of file

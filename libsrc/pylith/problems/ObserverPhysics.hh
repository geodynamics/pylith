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
 * @file libsrc/problems/ObserverPhysics.hh
 *
 * @brief Observer of physics (e.g., material, boundary condition, or interface condition).
 */

#if !defined(pylith_problems_observerphysics_hh)
#define pylith_problems_observerphysics_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/problems/Observer.hh" // ISA Observer

#include "pylith/feassemble/feassemblefwd.hh" // USES PhysicsImplementation
#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES PylithReal, PylithInt

class pylith::problems::ObserverPhysics : public pylith::problems::Observer {
    friend class TestObserverPhysics; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserverPhysics(void);

    /// Destructor
    virtual ~ObserverPhysics(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set physics implementation to observe.
     *
     * @param[in] physics Physics implementation to observe.
     */
    void setPhysicsImplementation(const pylith::feassemble::PhysicsImplementation* const physics);

    /** Set time scale.
     *
     * @param[in] value Time scale for dimensionalizing time.
     */
    virtual
    void setTimeScale(const PylithReal value) = 0;

    /** Verify observer is compatible with solution.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

    /** Receive update (subject of observer).
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     * @param[in] notification Type of notification.
     */
    virtual
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution,
                const pylith::problems::Observer::NotificationType notification) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////////////////////////////
protected:

    const pylith::feassemble::PhysicsImplementation* _physics; ///< Physics implementation to observe.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ObserverPhysics(const ObserverPhysics&); ///< Not implemented.
    const ObserverPhysics& operator=(const ObserverPhysics&); ///< Not implemented

}; // ObserverPhysics

#endif // pylith_problems_observerphysics_hh

// End of file

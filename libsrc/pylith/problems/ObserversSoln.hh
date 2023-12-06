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
 * @file libsrc/problems/ObserversSoln.hh
 *
 * @brief Collection of observers for an object.
 */

#if !defined(pylith_problems_observerssoln_hh)
#define pylith_problems_observerssoln_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/problems/Observer.hh" // USES NotificationType

#include "pylith/utils/GenericComponent.hh" // ISA GenericComponent

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES PylithReal, PylithInt

#include <set> // USES std::set

class pylith::problems::ObserversSoln : public pylith::utils::GenericComponent {
    friend class TestObserversSoln; // unit testing
    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserversSoln(void);

    /// Destructor
    virtual ~ObserversSoln(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Register observer to receive notifications.
     *
     * ObserversSoln are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(pylith::problems::ObserverSoln* observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(pylith::problems::ObserverSoln* observer);

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
     * @param[in] notification Type of notification.
     */
    void notifyObservers(const PylithReal t,
                         const PylithInt tindex,
                         const pylith::topology::Field& solution,
                         const pylith::problems::Observer::NotificationType notification);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Comparison function for keeping set of observers in order.
     *
     * @param[in] a Solution observer a.
     * @param[in] b Solution observer b.
     * @returns True if a->index < b->index, False otherwise.
     */
    struct _compare {
        bool operator()(const ObserverSoln* a,
                        const ObserverSoln* b) const;

    };

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    typedef std::set<pylith::problems::ObserverSoln*>::iterator iterator; ///< Iterator.
    std::set<pylith::problems::ObserverSoln*, _compare> _observers; ///< Subscribers of updates.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserversSoln(const ObserversSoln&); ///< Not implemented.
    const ObserversSoln& operator=(const ObserversSoln&); ///< Not implemented

};

// ObserversSoln

#endif // pylith_problems_observerssoln_hh

// End of file

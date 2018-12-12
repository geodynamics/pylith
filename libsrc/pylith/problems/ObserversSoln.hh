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
 * @file libsrc/problems/ObserversSoln.hh
 *
 * @brief Collection of observers for an object.
 */

#if !defined(pylith_problems_observerssoln_hh)
#define pylith_problems_observerssoln_hh

#include "problemsfwd.hh" // forward declarations

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
                         const pylith::topology::Field& solution);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    typedef std::set<pylith::problems::ObserverSoln*>::iterator iterator; ///< Iterator.
    std::set<pylith::problems::ObserverSoln*> _observers; ///< Subscribers of updates.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserversSoln(const ObserversSoln&); ///< Not implemented.
    const ObserversSoln& operator=(const ObserversSoln&); ///< Not implemented

};

// ObserversSoln

#endif // pylith_problems_observerssoln_hh

// End of file

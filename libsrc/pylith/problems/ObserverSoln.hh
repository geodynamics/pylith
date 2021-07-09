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
 * @file libsrc/problems/ObserverSoln.hh
 *
 * @brief ObserverSoln of subject.
 */

#if !defined(pylith_problems_observersoln_hh)
#define pylith_problems_observersoln_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES PylithReal, PylithInt

class pylith::problems::ObserverSoln {
    friend class TestObserverSoln; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor.
    ObserverSoln(void);

    /// Destructor
    virtual ~ObserverSoln(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

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
     * @param[in] infoOnly Flag is true if this update is before solution is available (e.g., after initialization).
     */
    virtual
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution) = 0;

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ObserverSoln(const ObserverSoln&); ///< Not implemented.
    const ObserverSoln& operator=(const ObserverSoln&); ///< Not implemented

}; // ObserverSoln

#endif // pylith_problems_observersoln_hh

// End of file

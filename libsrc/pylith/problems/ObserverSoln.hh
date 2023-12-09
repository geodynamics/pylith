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
 * @file libsrc/problems/ObserverSoln.hh
 *
 * @brief ObserverSoln of subject.
 */

#if !defined(pylith_problems_observersoln_hh)
#define pylith_problems_observersoln_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/problems/Observer.hh" // ISA Observer

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES PylithReal, PylithInt

class pylith::problems::ObserverSoln : public pylith::problems::Observer {
    friend ObserversSoln; ///< Access to ordering index.
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
     * @param[in] notification Type of notification.
     */
    virtual
    void update(const PylithReal t,
                const PylithInt tindex,
                const pylith::topology::Field& solution,
                const NotificationType notification) = 0;

    // PRIVATE ////////////////////////////////////////////////////////////////////////////////////
private:

    size_t index; ///< Index for keeing set of observers ordered.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    ObserverSoln(const ObserverSoln&); ///< Not implemented.
    const ObserverSoln& operator=(const ObserverSoln&); ///< Not implemented

}; // ObserverSoln

#endif // pylith_problems_observersoln_hh

// End of file

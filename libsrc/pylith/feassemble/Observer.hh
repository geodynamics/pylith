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
 * @file libsrc/feassemble/Observer.hh
 *
 * @brief Observer of subject.
 */

#if !defined(pylith_feassemble_observer_hh)
#define pylith_feassemble_observer_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/types.hh" // USES PylithReal, PylithInt

// Observer --------------------------------------------------------
/// Observer of subject. Receives updates of solution.
class pylith::feassemble::Observer {
    friend class TestObserver;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor.
    Observer(void);

    /// Destructor
    virtual ~Observer(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

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
                const pylith::topology::Field& solution,
                const bool infoOnly=false) = 0;

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    Observer(const Observer&);   ///< Not implemented.
    const Observer& operator=(const Observer&);   ///< Not implemented

}; // Observer

#endif // pylith_feassemble_observer_hh


// End of file

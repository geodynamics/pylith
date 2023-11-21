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
 * @file libsrc/problems/Observer.hh
 *
 * @brief Observer of subject.
 */

#if !defined(pylith_problems_observer_hh)
#define pylith_problems_observer_hh

#include "problemsfwd.hh" // forward declarations

class pylith::problems::Observer {
    friend class TestObserver; // unit testing

    // PUBLIC ENUMS ///////////////////////////////////////////////////////////////////////////////
public:

    enum NotificationType {
        DIAGNOSTIC=0,
        SOLUTION=1,
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Observer(void);

    /// Destructor
    virtual ~Observer(void);

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    Observer(const Observer&); ///< Not implemented.
    const Observer& operator=(const Observer&); ///< Not implemented

}; // Observer

#endif // pylith_problems_observer_hh

// End of file

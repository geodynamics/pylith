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

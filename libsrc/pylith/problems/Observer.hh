// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/problems/problemsfwd.hh" // forward declarations

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

// End of file
